import pandas as pd
from enum import Enum
from collections import Counter
from abc import abstractmethod
from glowingmeme.clients.clients import Clients
from protocols.protocol_7_2.reports import Program, Assembly


class ReportedOutcomeEnum(Enum):
    REPORTED = "reported"
    NOT_REPORTED = "not_reported"


class BuildDataset:
    DATASET_COLUMN_VALUES = ["id", "chromosome", "start", "end", "alt", "ref", "assembly", "case_id", "rs_id",
                             "age", "sex", "zigosity", "tier", "mode_of_inheritance", "consequence_type",
                             "MAF", "CADD_score", "type", "clinVar", "DisGeNET", "PhastCons", "phylop", "GERP",
                             "program", "mother_ethnic_origin", "father_ethnic_origin",
                             "gel_variant_acmg_classification", "case_solved_family", "phenotypes_solved",
                             "interpretation_message", "dict_extra_scores", "reported_outcome"]

    @abstractmethod
    def build_dataset(self):
        """
        This method should initialize the process of building a dataset and return a pandas DataFrame of it.
        :return:
        """
        pass


class BuildDatasetCVA(BuildDataset):

    def __init__(self):
        """
        Clients used to query services.
        :param cipapi:
        :param cva:
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
        reported_variant_list, non_reported_variant_list = self._query_cva_archived_positive_cases()
        self.main_dataset_df = self._create_dataframe_from_list(reported_variant_list + non_reported_variant_list)

    def _create_dataframe_from_list(self, list_to_add):
        """
        This method adds a given list to the object's main_dataset_df.
        :return:
        """
        dataset_df = pd.DataFrame(list_to_add, columns=self.DATASET_COLUMN_VALUES)
        return dataset_df

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
            case_id = "{identifier}-{version}".format(identifier=case.get("identifier", ""),
                                                      version=str(case.get("version", "")))
            sex = case.get("probandSex", None)
            program = case.get("program", None)
            tiered_variants = case.get("tieredVariants", {})
            age = case.get("probandEstimatedAgeAtAnalysis", None)
            classified_variants = case.get("classifiedVariants", {})
            interpretation_message = case.get("interpretation", None)

            for variant in case.get("reportedVariants", []):

                tier = self._get_variant_info(variant, tiered_variants)
                variant_acmg_classification = self._get_variant_info(variant, classified_variants)
                # variant corresponds to the queriable CVA id
                reported_variant_list.append([variant, None, None, None,  None, None, assembly,
                                              case_id, None, age, sex, None, tier,
                                              None, None, None, None, None, None,
                                              None, None, None, None, program, None, None, variant_acmg_classification,
                                              None, None, interpretation_message, None,
                                              ReportedOutcomeEnum.REPORTED.value])

            non_reported_variants = self._subtract_lists(case.get("reportedVariants", []),
                                                         case.get("allVariants", [])
                                                         )
            for variant in non_reported_variants:
                tier = self._get_variant_info(variant, tiered_variants)
                # variant corresponds to the queriable CVA id
                non_reported_variant_list.append([variant, None, None, None,  None, None, assembly,
                                              case_id, None, age, sex, None, tier,
                                              None, None, None, None, None, None,
                                              None, None, None, None, program, None, None, None,
                                              None, None, interpretation_message, None,
                                              ReportedOutcomeEnum.NOT_REPORTED.value])

        return reported_variant_list, non_reported_variant_list

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