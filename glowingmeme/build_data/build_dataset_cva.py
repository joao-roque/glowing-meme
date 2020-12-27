import os
from collections import Counter
from multiprocessing.dummy import Pool as ThreadPool

from protocols.protocol_7_2.reports import Program, Assembly

from glowingmeme.clients.clients import renew_access_token
from glowingmeme.build_data.build_dataset import BuildDataset
from glowingmeme.build_data.build_dataset import ReportedOutcomeEnum
from glowingmeme.build_data.variant_entry_info import VariantEntryInfo


class BuildDatasetCVA(BuildDataset):

    _CHROMOSOME = "chr"

    def __init__(self):
        """
        This is the first BuildDataset object to be called since it will fetch the relevant cases from CVA from which
        the remaining data will be fetched for.
        :return:
        """
        super().__init__()

    def build_dataset(self):
        """
        This method starts the process to build the Dataset Based on CVA queries.
        :return:
        """
        (
            reported_variant_list,
            non_reported_variant_list,
        ) = self._query_cva_archived_cases()
        self.main_dataset = reported_variant_list + non_reported_variant_list
        self._set_dataset_index_helper_by_attribute("id")
        self._fetch_specific_variant_information()

    @renew_access_token
    def _query_cva_archived_cases(self):
        """
        This method queries all the CVA cases that were archived with a positive result. It build
        :return: reported_variant_list, non_reported_variant_list
        """
        reported_variant_list = []
        non_reported_variant_list = []

        cases_iterator = self.cva_cases_client.get_cases(
            program=Program.rare_disease,
            assembly=Assembly.GRCh38,
            caseStatuses=["ARCHIVED_POSITIVE", "ARCHIVED_NEGATIVE"],
            include_all=False,
            max_results=1,
        )

        for case in cases_iterator:

            # since the variants belong to the same case, both the reported and non reported ones will have
            # some similar information, e.g. population
            assembly = case.get("assembly", None)
            case_id = "{identifier}-{version}".format(
                identifier=case.get("identifier", ""),
                version=str(case.get("version", "")),
            )
            sex = case.get("probandSex", None)
            program = case.get("program", None)
            tiered_variants = case.get("tieredVariants", {})
            age = case.get("probandEstimatedAgeAtAnalysis", None)
            classified_variants = case.get("classifiedVariants", {})
            interpretation_message = str(case.get("interpretation", None)).encode("utf-8")

            for variant in case.get("reportedVariants", []):
                tier = self._get_variant_info(variant, tiered_variants)
                variant_acmg_classification = self._get_variant_info(
                    variant, classified_variants
                )
                # variant corresponds to the queryable CVA id
                variant_info = VariantEntryInfo(
                    **{
                        "id": variant,
                        "assembly": assembly,
                        "case_id": case_id,
                        "age": age,
                        "sex": sex,
                        "tier": tier,
                        "program": program,
                        "gel_variant_acmg_classification": variant_acmg_classification,
                        "reported_outcome": ReportedOutcomeEnum.REPORTED.value,
                    }
                )

                reported_variant_list.append(variant_info)

            non_reported_variants = self._subtract_lists(
                case.get("reportedVariants", []), case.get("allVariants", [])
            )
            for variant in non_reported_variants:
                tier = self._get_variant_info(variant, tiered_variants)
                # variant corresponds to the queryable CVA id

                variant_info = VariantEntryInfo(
                    **{
                        "id": variant,
                        "assembly": assembly,
                        "case_id": case_id,
                        "age": age,
                        "sex": sex,
                        "tier": tier,
                        "program": program,
                        "interpretation_message": interpretation_message,
                        "reported_outcome": ReportedOutcomeEnum.NOT_REPORTED.value,
                    }
                )
                non_reported_variant_list.append(variant_info)

        return reported_variant_list, non_reported_variant_list

    def _fetch_specific_variant_information(self):
        """
        This method uses the CVA variant client and the previously fetched variants to provide more info for them.
        :return:
        """
        all_unique_variants = self.dataset_index_helper.keys()

        # threading available in client keeps breaking.
        # Implementing it here instead
        pool = ThreadPool(os.cpu_count())
        pool.daemon = True

        pool.map(self._query_and_fill_variant_object, all_unique_variants)

    @renew_access_token
    def _query_and_fill_variant_object(self, variant_id):
        """
        This method queries one variant ID, and updates the corresponding variant info object with the new info.
        :param variant_id:
        :return:
        """

        variant_wrapper = self.cva_variants_client.get_variant_by_id(variant_id)

        for variant in variant_wrapper.variants:

            if variant.assembly == self._ASSEMBLY_38 and variant.annotation:
                for variant_info_object in self.dataset_index_helper[
                    variant_wrapper.id
                ]:

                    variant_type = variant.smallVariantType
                    if variant.variantType:
                        variant_type = variant.variantType

                    # rebuilding list with new values fetched
                    variant_info_object.update_object(
                        **{
                            "chromosome": self._CHROMOSOME + variant.annotation.chromosome,
                            "start": variant.annotation.start,
                            "end": variant.annotation.start + len(variant.annotation.reference),
                            "alt": variant.annotation.alternate,
                            "ref": variant.annotation.reference,
                            "rs_id": variant.annotation.id,
                            "consequence_type": self._get_sequence_ontology_terms(
                                variant.annotation.consequenceTypes
                            ),
                            "biotypes": self._get_biotypes(
                                variant.annotation.consequenceTypes
                            ),
                            "population_frequency": self._get_population_frequency(
                                variant_info_object.sex,
                                variant.annotation.populationFrequencies,
                            ),
                            "type": variant_type,
                            "PhastCons": self._get_conservation_score_from_source(
                                self._PHAST_CONS, variant.annotation.conservation,
                            ),
                            "phylop": self._get_conservation_score_from_source(
                                self._PHYLOP, variant.annotation.conservation
                            ),
                        }
                    )

    def _get_population_frequency(self, sex, population_frequencies_list):
        """
        This method, from a given PopulationsFrequencies object and sex of participant, extracts the most relevant
        population frequency. For now it only uses sex for its logic.
        :return:
        """
        population_all_frequency = None
        population_sex_frequency = None
        if population_frequencies_list:
            for population_frequency in population_frequencies_list:
                if population_frequency.study == self._GNOMAD_GENOMES:
                    if population_frequency.population == sex:
                        population_sex_frequency = population_frequency.altAlleleFreq

                    if population_frequency.population == self._ALL:
                        population_all_frequency = population_frequency.altAlleleFreq

        if population_sex_frequency:
            return population_sex_frequency
        return population_all_frequency

    @staticmethod
    def _get_sequence_ontology_terms(consequence_types_list):
        """
        This method returns a list of sequence ontology terms given a list of consequenceType objects.
        :return:
        """
        sequence_names = []
        for consequence_type in consequence_types_list:
            for sequence_ontology in consequence_type.sequenceOntologyTerms:
                if sequence_ontology.name:
                    sequence_names.append(sequence_ontology.name)

        if sequence_names:
            return ",".join(sequence_names)
        return None

    @staticmethod
    def _get_biotypes(consequence_types_list):
        """
        This method returns a list of biotypes given a list of consequenceType objects.
        :return:
        """
        biotypes = [
            consequence_type.biotype
            for consequence_type in consequence_types_list
            if consequence_type.biotype
        ]

        if biotypes:
            return ",".join(biotypes)
        return None

    @staticmethod
    def _get_conservation_score_from_source(source, conservation):
        """
        This method returns the required score given a conservation object list.
        :param conservation:
        :return:
        """
        if conservation:
            for conservation_score in conservation:
                if conservation_score.source == source:
                    return conservation_score.score
        return None

    @staticmethod
    def _get_variant_info(variant, info_dict):
        """
        This method returns the tier for a given variant given a dictionary of the variant tiers.
        :param variant:
        :param info_dict:
        :return:
        """
        info = None
        for key in info_dict.keys():
            if variant in set(info_dict[key]):
                info = key
        return info

    @staticmethod
    def _subtract_lists(array_1, array_2):
        """
        This method subtracts elements of array_2 from array_1
        :return:
        """
        subtracted_list = set((Counter(array_2) - Counter(array_1)).elements())
        return subtracted_list

