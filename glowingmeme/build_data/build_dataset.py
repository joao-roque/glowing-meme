import pandas as pd
from enum import Enum
from collections import Counter
from abc import abstractmethod
from glowingmeme.clients.clients import Clients
from protocols.protocol_7_2.reports import Program, Assembly


class ReportedOutcomeEnum(Enum):
    REPORTED = "reported"
    NOT_REPORTED = "not_reported"


class BuildDataset(object):

    @abstractmethod
    def build_dataset(self):
        """
        This method should initialize the process of building a dataset and return a pandas DataFrame of it.
        :return:
        """
        pass


class BuildDatasetCVA(BuildDataset):
    DATASET_COLUMN_VALUES = ["id", "chromosome", "start", "end", "build", "assembly", "case_id", "rs_id", "population",
                             "age", "sex", "zigosity", "biotype", "tier", "mode_of_inheritance", "consequence_type",
                             "MAF", "CADD_score", "type", "clinVar", "DisGeNET", "PhastCons", "phylop", "GERP",
                             "gene_name", "program", "reported_outcome"]

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
        self.main_dataset_df = pd.DataFrame(columns=self.DATASET_COLUMN_VALUES)

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
            age = case.get("probandEstimatedAgeAtAnalysis", None)
            tiered_variants = case.get("tieredVariants", {})

            for variant in case.get("reportedVariants", []):

                tier = self._get_variant_tier(variant, tiered_variants)
                # variant corresponds to the queriable CVA id
                reported_variant_list.append([variant, None, None, None, assembly,
                                              case_id,
                                              None, age, sex, None, None, tier,
                                              None, None, None, None, None, None, None,
                                              None, None, None, None, program,
                                              ReportedOutcomeEnum.REPORTED.value])

            non_reported_variants = self._subtract_lists(case.get("reportedVariants", []),
                                                         case.get("allVariants", [])
                                                         )
            for variant in non_reported_variants:

                tier = self._get_variant_tier(variant, tiered_variants)
                # variant corresponds to the queriable CVA id
                non_reported_variant_list.append([variant, None, None, None, assembly,
                                                  case_id,
                                                  None, age, sex, None, None, tier,
                                                  None, None, None, None, None, None, None,
                                                  None, None, None, None, program,
                                                  ReportedOutcomeEnum.NOT_REPORTED.value])

        return reported_variant_list, non_reported_variant_list

    def _get_variant_tier(self, variant, tiered_variants):
        """
        This method returns the tier for a given variant given a dictionary of the variant tiers.
        :param variant:
        :param tiered_variants:
        :return:
        """
        tier = None
        for key in tiered_variants.keys():
            if variant in set(tiered_variants[key]):
                tier = key
        return tier

    @staticmethod
    def _subtract_lists(array_1, array_2):
        """
        This method subtracts elements of array_2 from array_1
        :return:
        """
        subtracted_list = set((Counter(array_2) - Counter(array_1)).elements())
        return subtracted_list

    def get_variants_used_in_diagnosis(self):
        """
        This methods queries Genomics England services
        to gather variants that were used in the diagnosis by the clinician.
        :return:
        """
        pass
