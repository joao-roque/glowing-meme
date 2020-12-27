from multiprocessing.dummy import Pool as ThreadPool

from glowingmeme.clients.clients import renew_access_token
from glowingmeme.build_data.build_dataset import BuildDataset


class BuildDatasetCellbase(BuildDataset):

    _POST = "post"
    _RESULT_FIELD = "result"
    _PHAST_CONS = "phastCons"
    _VARIANT_ID_SEPARATOR = ":"
    _CELLBASE_QUERY_BATCH_SIZE = 1000
    _VARIANT_TRAIT_ASSOCIATION = "variantTraitAssociation"

    _INCLUDE_LIST = [
        _VARIANT_TRAIT_ASSOCIATION,
        "gwas",
        "functionalScore"
    ]

    def __init__(self, cipapi_built_dataset):
        """
        This method takes in its precursor dataset, which is the one updated by the Cipapi Dataset builder.
        :param cipapi_built_dataset:
        """
        super().__init__()
        self.main_dataset = cipapi_built_dataset
        self.dataset_index_helper = None

    def build_dataset(self):
        """
        This method starts updating the variant entries with Cellbase info.
        :return:
        """
        self._populate_cellbase_query_id()
        self._set_dataset_index_helper_by_attribute("cellbase_query_id")
        self._annotate_variants_in_batches()

    def _build_variant_id(self, variant_info_object):
        """
        This method creates a variants ids to query Cellbase from the variant object
        :param variant_info_object:
        :return:
        """
        # for some reason, sometimes CNVs are included with the <CN_> notation. We skip those here
        if "<" in variant_info_object.ref or "<" in variant_info_object.alt:
            return ""

        return self._VARIANT_ID_SEPARATOR.join([variant_info_object.chromosome,
                                                str(variant_info_object.start),
                                                variant_info_object.ref,
                                                variant_info_object.alt])

    def _populate_cellbase_query_id(self):
        """
        This method populates the cellbase_query_id for each Variant Info object.
        :return:
        """
        for variant_info in self.main_dataset:
            variant_info.update_object(**{"cellbase_query_id": self._build_variant_id(variant_info)})

    def _annotate_variants_in_batches(self):
        """
        We build batches of variants to query Cellbase with, so that we don't send too big of a request which
        threatens to send the service down.
        :return:
        """

        # TODO this needs to be changed to divide by unique variants to be queried
        list_of_batches = []
        variant_ids_to_query = [variant_id for variant_id in self.dataset_index_helper.keys() if variant_id]
        for list_chunk in range(0, len(variant_ids_to_query), self._CELLBASE_QUERY_BATCH_SIZE):
            list_of_batches.append(variant_ids_to_query[list_chunk:list_chunk + self._CELLBASE_QUERY_BATCH_SIZE])

        # we are putting a hard cap of threads here to not overload Cellbase
        pool = ThreadPool(8)
        pool.map(self._call_cellbase_with_batch, list_of_batches)

    @renew_access_token
    def _call_cellbase_with_batch(self, variant_ids_to_query):
        """
        This method takes a list of VariantInfo objects. It will query Cellbase and update the object with the new info.
        :param variants_to_query:
        :return:
        """

        for variant_id_list_index, response in enumerate(self.cellbase_client.get_annotation(variant_ids_to_query,
                                                                                  method=self._POST,
                                                                                  # include=self._INCLUDE_LIST,
                                                                                  assembly=self._ASSEMBLY_38,
                                                                                  )):
            if response[self._RESULT_FIELD]:
                cellbase_variant_info = response[self._RESULT_FIELD][0]

                # now we fill in all the variant information for these variants
                for variant_info in self.dataset_index_helper[variant_ids_to_query[variant_id_list_index]]:

                    # There should always be only one response. However if there are more here, we ignore them
                    variant_info.update_object(**{
                        "CADD_score": cellbase_variant_info,
                        "clinVar": self._get_clinvar_classification(cellbase_variant_info[self._VARIANT_TRAIT_ASSOCIATION]),

                    })

    def _get_clinvar_classification(self, variant_trait_association):
        """
        From the variantTraitAssociation response from Cellbase, extract ClinVar information
        :param variant_trait_association:
        :return:
        """

        if "clinvar" in variant_trait_association:
            # if there are only numbers in the accession, that is the one we want, as opposed to the individual records.
            for record in variant_trait_association["clinvar"]:
                if record["accession"].isdigit():
                    return record["clinicalSignificance"]
        return
