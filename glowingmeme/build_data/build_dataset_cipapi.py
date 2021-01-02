import os
from operator import attrgetter
from multiprocessing.dummy import Pool as ThreadPool

from glowingmeme.clients.clients import renew_access_token
from glowingmeme.build_data.build_dataset import BuildDataset


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
        super().__init__()
        self.main_dataset = cva_built_dataset
        self.dataset_index_helper = None

    def build_dataset(self):
        """
        Start building the dataset.
        :return:
        """
        self._fetch_cipapi_data()

    def _fetch_cipapi_data(self):
        """
        This method queries cipapi for more data.
        :return:
        """

        self._set_dataset_index_helper_by_attribute("case_id")
        case_id_list = self.dataset_index_helper.keys()

        pool = ThreadPool(os.cpu_count())
        pool.daemon = True

        pool.map(self._query_data_for_case, case_id_list)

    @renew_access_token
    def _query_data_for_case(self, case_id):
        """
        This method queries all the data for a given cipapi case.
        :param case_id:
        :return:
        """

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

        fast_lookup_dict = self._create_fast_lookup_dict(
            genomics_england_interpreted_genome
        )

        # we can now fill in every variant for this case with the relevant information
        for variant_entry_info in self.dataset_index_helper[case_id]:

            fast_lookup_key_name = "{chr}_{start}".format(
                chr=variant_entry_info.chromosome.replace("chr", ""),
                start=variant_entry_info.start,
            )

            # get corresponding variant from interpreted genome
            variant_in_genome = None
            if fast_lookup_key_name in fast_lookup_dict:
                variant_in_genome = fast_lookup_dict[fast_lookup_key_name][0]

            if variant_in_genome:
                (
                    proband_zygosity,
                    mother_zygosity,
                    father_zygosity,
                ) = self._get_family_variant_zygosity(
                    pedigree=interpretation_request.pedigree, variant=variant_in_genome,
                )

                variant_entry_info.zygosity_proband = proband_zygosity
                variant_entry_info.zygosity_mother = mother_zygosity
                variant_entry_info.zygosity_father = father_zygosity

                variant_entry_info.mode_of_inheritance = variant_in_genome.reportEvents[
                    0
                ].modeOfInheritance
                variant_entry_info.segregation_pattern = variant_in_genome.reportEvents[
                    0
                ].segregationPattern
                variant_entry_info.penetrance = variant_in_genome.reportEvents[
                    0
                ].penetrance

                variant_entry_info.mother_ethnic_origin = (
                    proband.ancestries.mothersEthnicOrigin
                )
                variant_entry_info.father_ethnic_origin = (
                    proband.ancestries.fathersEthnicOrigin
                )

                if (
                    latest_report.exit_questionnaire
                    and len(
                        latest_report.exit_questionnaire.exit_questionnaire_data[
                            "variantGroupLevelQuestions"
                        ]
                    )
                    >= 1
                ):

                    variant_entry_info.case_solved_family = latest_report.exit_questionnaire.exit_questionnaire_data[
                        "familyLevelQuestions"
                    ][
                        "caseSolvedFamily"
                    ]

                    variant_entry_info.phenotypes_solved = latest_report.exit_questionnaire.exit_questionnaire_data[
                        "variantGroupLevelQuestions"
                    ][
                        -1
                    ][
                        "phenotypesSolved"
                    ]

                    variant_entry_info.actionability = latest_report.exit_questionnaire.exit_questionnaire_data[
                        "variantGroupLevelQuestions"
                    ][
                        -1
                    ][
                        "actionability"
                    ]

    def _create_fast_lookup_dict(self, interpreted_genome):
        """
        This method takes all the small variant objects from an interpreted genome and puts them in a dictionary
        index by chr_start. This makes it extremely quicker to query for variants later, since we can reduce the
        searchable list by that index.
        :param interpreted_genome:
        :return:
        """

        fast_lookup_variant_dict = {}
        for variant in interpreted_genome.interpretation_request_payload.variants:
            key_name = "{chr}_{start}".format(
                chr=variant.variantCoordinates.chromosome,
                start=variant.variantCoordinates.position,
            )

            if key_name in fast_lookup_variant_dict:
                fast_lookup_variant_dict[key_name].append(variant)
            else:
                fast_lookup_variant_dict[key_name] = [variant]

        return fast_lookup_variant_dict

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
