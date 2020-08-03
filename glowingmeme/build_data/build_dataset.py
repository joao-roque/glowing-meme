import pandas as pd
from glowingmeme.clients.clients import Clients
from protocols.protocol_7_2.reports import Program, Assembly


class BuildDataset:

    DATASET_COLUMN_VALUES = ["id", "chromosome", "start", "end", "build", "case_id", "rs_id", "population", "age",
                             "sex", "zigosity", "biotype", "tier", "mode_of_inheritance", "consequence_type", "MAF",
                             "CADD_score", "type", "clinVar", "DisGeNET", "PhastCons", "phylop", "GERP", "gene_name",
                             "program", "reported_outcome"]

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
        This method queries all the CVA cases that were archived with a positive result.
        :return:
        """

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
            program = case.get("program", None)
            age = case.get("probandEstimatedAgeAtAnalysis", None)
            sex = case.get("probandSex", None)

            for variant in case.get("reportedVariants", []):

                #self.main_dataset_df.iloc[len(self.main_dataset_df.index) + 1] =
                pass



    def get_variants_used_in_diagnosis(self):
        """
        This methods queries Genomics England services
        to gather variants that were used in the diagnosis by the clinician.
        :return:
        """
        pass
