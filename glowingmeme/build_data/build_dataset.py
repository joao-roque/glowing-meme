import csv
from enum import Enum
from abc import abstractmethod

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

        # predefining clients that will be used
        self.cva_client = None
        self.cipapi_client = None
        self.cellbase_client = None
        self.cva_cases_client = None
        self.cva_variants_client = None

        self.start_clients()

        # IMPORTANT DESCRIPTION OF DATASET
        # The dataset IS composed of variants that are associated with a specific case. This means that the same
        # variant can appear multiple times along the dataset as long as it does not contain repeated information
        # and outcomes.

        # this is a list of VariantEntryInfo objects
        self.main_dataset = []

        # this helper can be redefined by which attribute need by calling _set_dataset_index_helper_by_attribute
        # NOTE: This dictionary is a reference to the original objects in the main_dataset list. If a change is made
        # in these object, the original objects in the main_dataset will also change.
        self.dataset_index_helper = {}

    @abstractmethod
    def build_dataset(self):
        """
        This method should initialize the process of building a dataset.
        :return:
        """
        pass

    def start_clients(self):
        """
        Start all the required clients.
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
                new_key = getattr(variant_info_object, dataset_key)
                if new_key in dataset_by_key:
                    dataset_by_key[new_key].append(variant_info_object)
                else:
                    dataset_by_key[new_key] = [variant_info_object]
        self.dataset_index_helper = dataset_by_key

    def save_data_to_csv(self, file_name):
        """
        This method takes the main dataset that was created and saves it to a csv
        :param file_name:
        :return:
        """
        with open(file_name, "w") as variant_entries_file:

            variant_entries_csv = csv.writer(variant_entries_file, delimiter=",")

            # we start by adding the header, which is always the variant_info_values from the object
            variant_entries_csv.writerow(VariantEntryInfo.VARIANT_INFO_VALUES)

            for variant_entry in self.main_dataset:
                variant_entries_csv.writerow(list(variant_entry))
