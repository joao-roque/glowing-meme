import pandas as pd
from enum import Enum
from collections import Counter
from abc import abstractmethod
from glowingmeme.clients.clients import Clients
from multiprocessing.dummy import Pool as ThreadPool
from protocols.protocol_7_2.reports import Program, Assembly

class ReportedOutcomeEnum(Enum):
    REPORTED = "reported"
    NOT_REPORTED = "not_reported"


class BuildDataset:
    _ALL = "ALL"
    _MALE = "MALE"
    _FEMALE = "FEMALE"
    _PHYLOP = "phylop"
    _ASSEMBLY_38 = "GRCh38"
    _PHAST_CONS = "phastCons"
    _GNOMAD_GENOMES = 'GNOMAD_GENOMES'

    _DATASET_COLUMN_VALUES = [
        "id",
        "chromosome",
        "start",
        "end",
        "alt",
        "ref",
        "assembly",
        "case_id",
        "rs_id",
        "age",
        "sex",
        "zigosity",
        "tier",
        "mode_of_inheritance",
        "consequence_type",
        "biotypes",
        "population_frequency",
        "CADD_score",
        "type",
        "clinVar",
        "DisGeNET",
        "PhastCons",
        "phylop",
        "GERP",
        "program",
        "mother_ethnic_origin",
        "father_ethnic_origin",
        "gel_variant_acmg_classification",
        "case_solved_family",
        "phenotypes_solved",
        "interpretation_message",
        "dict_extra_scores",
        "reported_outcome",
    ]

    @abstractmethod
    def build_dataset(self):
        """
        This method should initialize the process of building a dataset and return a pandas DataFrame of it.
        :return:
        """
        pass

    def _create_dataset_list(self, **dataset_column_values):
        """
        This method creates a list with parameters from DATASET_COLUMN_VALUES from a given dict. If one of the
         parameters does not exist it will be set to None.
        :param dataset_column_values:
        :return:
        """
        dataset_list = []
        for key in self._DATASET_COLUMN_VALUES:
            if key in dataset_column_values:
                dataset_list.append(dataset_column_values[key])
            else:
                dataset_list.append(None)
        return dataset_list

    def _start_clients(self):
        """
        Auto restart of clients so tokens don't expire.
        :return:
        """
        (
            self.cipapi_client,
            self.cellbase_client,
            self.cva_client,
        ) = Clients().get_all_clients()
        self.cva_cases_client = self.cva_client.cases()
        self.cva_variants_client = self.cva_client.variants()


class BuildDatasetCVA(BuildDataset):
    def __init__(self):
        """
        Clients used to query services.
        """
        self._start_clients()

        # IMPORTANT DESCRIPTION OF DATASET
        # The dataset IS composed of variants that are associated with a specific case. This means that the same
        # variant can appear multiple times along the dataset as long as it does not contain repeated information
        # and outcomes.

        # initialized in build dataset
        self.main_dataset_df = None

    def build_dataset(self):
        """
        This method starts the process to build the Dataset Based on CVA queries.
        :return:
        """
        (
            reported_variant_list,
            non_reported_variant_list,
        ) = self._query_cva_archived_positive_cases()
        self.main_dataset_df = self._create_dataframe_from_list(
            reported_variant_list + non_reported_variant_list
        )
        self._fetch_specific_variant_information()


    def _create_dataframe_from_list(self, list_to_add):
        """
        This method adds a given list to the object's main_dataset_df.
        :return:
        """
        dataset_df = pd.DataFrame(list_to_add, columns=self._DATASET_COLUMN_VALUES)
        return dataset_df

    def _query_cva_archived_positive_cases(self):
        """
        This method queries all the CVA cases that were archived with a positive result. It build
        :return: reported_variant_list, non_reported_variant_list
        """
        reported_variant_list = []
        non_reported_variant_list = []

        cases_iterator = self.cva_cases_client.get_cases(
            program=Program.rare_disease,
            assembly=Assembly.GRCh38,
            caseStatuses="ARCHIVED_POSITIVE",
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
            interpretation_message = case.get("interpretation", None)

            for variant in case.get("reportedVariants", []):
                tier = self._get_variant_info(variant, tiered_variants)
                variant_acmg_classification = self._get_variant_info(
                    variant, classified_variants
                )
                # variant corresponds to the queryable CVA id
                variant_info = self._create_dataset_list(
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

                variant_info = self._create_dataset_list(
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
        all_unique_variants = set(self.main_dataset_df["id"].tolist())

        # threading available in client keeps breaking.
        # Implementing it here instead
        pool = ThreadPool(8)
        all_unique_variants_fetched = pool.map(self.cva_variants_client.get_variant_by_id, all_unique_variants)

        for variant_wrapper in all_unique_variants_fetched:
            if variant_wrapper.variants:
                for variant in variant_wrapper.variants:
                    if variant.assembly == self._ASSEMBLY_38 and variant.annotation:
                        variant_in_dataset = self.main_dataset_df.loc[self.main_dataset_df["id"] == variant_wrapper.id]

                        variant_type = variant.smallVariantType
                        if variant.variantType:
                            variant_type = variant.variantType

                        # rebuilding list with new values fetched
                        variant_info = self._create_dataset_list(
                            **{
                                "id": variant_in_dataset["id"].values[0],
                                "chromosome": variant.annotation.chromosome,
                                "start": variant.annotation.start,
                                "end": variant.annotation.end,
                                "alt": variant.annotation.alternate,
                                "ref": variant.annotation.reference,
                                "rs_id": variant.annotation.id,
                                "consequence_type": self._get_sequence_ontology_terms(
                                    variant.annotation.consequenceTypes),
                                "biotypes": self._get_biotypes(variant.annotation.consequenceTypes),
                                "population_frequency": self._get_population_frequency(
                                    variant_in_dataset["sex"].values[0], variant.annotation.populationFrequencies),
                                "type": variant_type,
                                "PhastCons": self._get_conservation_score_from_source(self._PHAST_CONS,
                                                                                      variant.annotation.conservation),
                                "phylop": self._get_conservation_score_from_source(self._PHYLOP,
                                                                                   variant.annotation.conservation),
                                "assembly": variant_in_dataset["assembly"].values[0],
                                "case_id": variant_in_dataset["case_id"].values[0],
                                "age": variant_in_dataset["age"].values[0],
                                "sex": variant_in_dataset["sex"].values[0],
                                "tier": variant_in_dataset["tier"].values[0],
                                "program": variant_in_dataset["program"].values[0],
                                "gel_variant_acmg_classification":
                                    variant_in_dataset["gel_variant_acmg_classification"].values[0],
                                "interpretation_message": variant_in_dataset["interpretation_message"].values[0],
                                "reported_outcome": variant_in_dataset["reported_outcome"].values[0],
                            }
                        )

                        self.main_dataset_df.loc[self.main_dataset_df["id"] == variant_wrapper.id] = variant_info

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

    def _get_sequence_ontology_terms(self, consequence_types_list):
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

    def _get_biotypes(self, consequence_types_list):
        """
        This method returns a list of biotypes given a list of consequenceType objects.
        :return:
        """
        biotypes = [consequence_type.biotype for consequence_type in consequence_types_list if consequence_type.biotype]

        if biotypes:
            return ",".join(biotypes)
        return None

    def _get_conservation_score_from_source(self, source, conservation):
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

    def _get_variant_info(self, variant, info_dict):
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


class BuildDatasetCipapi(BuildDataset):
    def __init__(self, cva_built_dataset):
        """
        This class takes as precursor a Pandas Dataframe with the columns defined in the parent class, in the variable
        DATASET_COLUMN_VALUES. It requires at least the columns case_id, assembly and variant details
        ("chromosome", "start", "end") to be populated, otherwise it won't be able to find this information in cipapi.
        :param cva_built_dataset:
        """
        pass

    def build_dataset(self):
        """
        Start building the dataset.
        :return:
        """
        pass
