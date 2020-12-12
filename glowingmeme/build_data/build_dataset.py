from enum import Enum
from abc import abstractmethod
from collections import Counter
from operator import attrgetter
from multiprocessing.dummy import Pool as ThreadPool

from protocols.protocol_7_2.reports import Program, Assembly

from glowingmeme.clients.clients import Clients
from glowingmeme.build_data.variant_entry_info import VariantEntryInfo


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
    _GNOMAD_GENOMES = "GNOMAD_GENOMES"

    def __init__(self):
        self._start_clients()
        # IMPORTANT DESCRIPTION OF DATASET
        # The dataset IS composed of variants that are associated with a specific case. This means that the same
        # variant can appear multiple times along the dataset as long as it does not contain repeated information
        # and outcomes.

        # this is a set of VariantEntryInfo objects. It is a set because we do not want repeated variants
        self.main_dataset = set()

        # this helper can be redefined by which attribute need by calling _set_dataset_index_helper_by_attribute
        # NOTE: This dictionary is a reference to the original objects in the main_dataset list. If a change is made
        # in these object, the original objects in the main_dataset will also change.
        self.dataset_index_helper = {}

    @abstractmethod
    def build_dataset(self):
        """
        This method should initialize the process of building a dataset and return a pandas DataFrame of it.
        :return:
        """
        pass

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

    def _set_dataset_index_helper_by_attribute(self, dataset_key):
        """
        In order to massively speed up querying specific VariantInfo objects of the main dataset, we here create
        a dictionary that will index said objects by a given key e.g variantId

        This dictionary is a reference to the original objects in the main_dataset list
        :return:
        """
        dataset_by_key = {}
        if dataset_key in VariantEntryInfo.VARIANT_INFO_VALUES:
            for variant_info_object in self.main_dataset:
                new_key = vars(variant_info_object)[dataset_key]
                if new_key in dataset_by_key:
                    dataset_by_key[new_key].append(variant_info_object)
                else:
                    dataset_by_key[new_key] = [variant_info_object]
        self.dataset_index_helper = dataset_by_key


class BuildDatasetCVA(BuildDataset):
    def __init__(self):
        """
        Clients used to query services.
        """
        BuildDataset.__init__(self)

    def build_dataset(self):
        """
        This method starts the process to build the Dataset Based on CVA queries.
        :return:
        """
        (
            reported_variant_list,
            non_reported_variant_list,
        ) = self._query_cva_archived_cases()
        self.main_dataset.update(reported_variant_list + non_reported_variant_list)
        self._set_dataset_index_helper_by_attribute("id")
        self._fetch_specific_variant_information()
        return self.main_dataset

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
            interpretation_message = case.get("interpretation", None)

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
        pool = ThreadPool(8)
        all_unique_variants_fetched = pool.map(
            self.cva_variants_client.get_variant_by_id, all_unique_variants
        )

        for variant_wrapper in all_unique_variants_fetched:
            if variant_wrapper.variants:
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
                                    "chromosome": variant.annotation.chromosome,
                                    "start": variant.annotation.start,
                                    "end": variant.annotation.end,
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
                                        self._PHAST_CONS,
                                        variant.annotation.conservation,
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


class BuildDatasetCipapi(BuildDataset):

    FATHER = "Father"
    MOTHER = "Mother"
    GENOMICS_ENGLAND_TIERING = "genomics_england_tiering"

    def __init__(self, cva_built_dataset):
        """
        This class takes as precursor a Pandas Dataframe with the columns defined in the parent class, in the variable
        DATASET_COLUMN_VALUES. It requires at least the columns case_id, assembly and variant details
        ("chromosome", "start", "end") to be populated, otherwise it won't be able to find this information in cipapi.
        :param cva_built_dataset:
        """
        BuildDataset.__init__(self)
        self.main_dataset = cva_built_dataset
        self.dataset_index_helper = None

    def build_dataset(self):
        """
        Start building the dataset.
        :return:
        """
        self._set_dataset_index_helper_by_attribute("case_id")
        self._fetch_cipapi_data()
        return self.main_dataset

    def _fetch_cipapi_data(self):
        """
        This method queries cipapi to complete some of the dataframe's fields.
        :return:
        """

        for case_id in self.dataset_index_helper.keys():
            case = case_id.split("-")[0]
            version = case_id.split("-")[1]
            interpretation_request = self.cipapi_client.get_case(
                case_id=case, case_version=version
            )

            # always pull the info from the latest report
            latest_report = self._get_latest_report(
                interpretation_request=interpretation_request
            )
            proband = [
                member
                for member in interpretation_request.pedigree.members
                if member.isProband
            ][0]
            genomics_england_interpreted_genome = self._get_gel_interpreted_genome(
                interpretation_request=interpretation_request
            )

            # we can now fill in every variant for this case with the relevant information
            for variant_entry_info in self.dataset_index_helper[case_id]:

                # get corresponding variant from interpreted genome
                variant_in_genome = self._get_variant_from_interpreted_genome(
                    interpreted_genome=genomics_england_interpreted_genome,
                    variant_entry_info=variant_entry_info,
                )

                if variant_in_genome:
                    (
                        proband_zygosity,
                        mother_zygosity,
                        father_zygosity,
                    ) = self._get_family_variant_zygosity(
                        pedigree=interpretation_request.pedigree,
                        variant=variant_in_genome,
                    )

                    variant_entry_info.zygosity_proband = proband_zygosity
                    variant_entry_info.zygosity_mother = mother_zygosity
                    variant_entry_info.zygosity_father = father_zygosity

                    variant_entry_info.mode_of_inheritance = (
                        variant_in_genome.reportEvents[0].modeOfInheritance
                    )
                    variant_entry_info.segregation_pattern = (
                        variant_in_genome.reportEvents[0].segregationPattern
                    )
                    variant_entry_info.penetrance = variant_in_genome.reportEvents[
                        0
                    ].penetrance

                    variant_entry_info.mother_ethnic_origin = (
                        proband.ancestries.mothersEthnicOrigin
                    )
                    variant_entry_info.father_ethnic_origin = (
                        proband.ancestries.fathersEthnicOrigin
                    )

                    variant_entry_info.case_solved_family = (
                        latest_report.exit_questionnaire.exit_questionnaire_data[
                            "familyLevelQuestions"
                        ]["caseSolvedFamily"]
                    )
                    variant_entry_info.phenotypes_solved = (
                        latest_report.exit_questionnaire.exit_questionnaire_data[
                            "variantGroupLevelQuestions"
                        ][-1]["phenotypesSolved"]
                    )
                    variant_entry_info.actionability = (
                        latest_report.exit_questionnaire.exit_questionnaire_data[
                            "variantGroupLevelQuestions"
                        ][-1]["actionability"]
                    )

    def _get_family_variant_zygosity(self, pedigree, variant):
        """
        Given a pedigree and a variant, this method will return the zygosity values for the family.
        :param pedigree:
        :param variant:
        :return: proband_zygosity, mother_zygosity, father_zygosity
        """

        mother_zygosity = None
        father_zygosity = None
        proband_zygosity = None

        (
            proband_participant_id,
            mother_participant_id,
            father_participant_id,
        ) = self._get_family_ids(pedigree=pedigree)

        for variant_call in variant.variantCalls:
            if variant_call.participantId == proband_participant_id:
                proband_zygosity = variant_call.zygosity
            elif variant_call.participantId == mother_participant_id:
                mother_zygosity = variant_call.zygosity
            elif variant_call.participantId == father_participant_id:
                father_zygosity = variant_call.zygosity

        return proband_zygosity, mother_zygosity, father_zygosity

    def _get_variant_from_interpreted_genome(
        self, interpreted_genome, variant_entry_info
    ):
        """
        This method will find the variant in the variant_entry_info in the given interpreted genome and return it.
        :param interpreted_genome:
        :param variant_entry_info:
        :return:
        """
        interpreted_genome_variant_list = [
            variant
            for variant in interpreted_genome.interpretation_request_payload.variants
            if variant.variantCoordinates.alternate == variant_entry_info.alt
            and variant.variantCoordinates.chromosome == variant_entry_info.chromosome
            and variant.variantCoordinates.position == variant_entry_info.start
            and variant.variantCoordinates.reference == variant_entry_info.ref
        ]

        if len(interpreted_genome_variant_list) >= 1:
            return interpreted_genome_variant_list[0]

        return None

    def _get_gel_interpreted_genome(self, interpretation_request):
        """
        Given an Interpretation Request, the interpreted genome that was created by genomics england tiering services
        will be returned.
        :param interpretation_request:
        :return:
        """
        return max(
            [
                interpreted_genome
                for interpreted_genome in interpretation_request.interpreted_genome
                if interpreted_genome.interpretation_request_payload.interpretationService
                == self.GENOMICS_ENGLAND_TIERING
            ],
            key=attrgetter("created_at"),
        )

    def _get_family_ids(self, pedigree):
        """
        Given a Pedigree object, this method will return the proband, father and mother participant ids.
        :param pedigree:
        :return: proband, mother, father participantIds
        """
        mother_participant_id = None
        father_participant_id = None
        proband_participant_id = None

        for member in pedigree.members:
            if member.isProband:
                proband_participant_id = member.participantId
            elif member.additionalInformation["relation_to_proband"] == self.FATHER:
                father_participant_id = member.participantId
            elif member.additionalInformation["relation_to_proband"] == self.MOTHER:
                mother_participant_id = member.participantId

        return proband_participant_id, mother_participant_id, father_participant_id

    @staticmethod
    def _get_latest_report(interpretation_request):
        """
        From a given interpretation request, the latest clinical report is pulled
        :param interpretation_request:
        :return:
        """

        return max(
            [
                clinical_report
                for clinical_report in interpretation_request.clinical_report
            ],
            key=attrgetter("created_at"),
        )
