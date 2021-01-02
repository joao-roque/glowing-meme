import os
from multiprocessing.dummy import Pool as ThreadPool

from glowingmeme.clients.clients import renew_access_token
from glowingmeme.build_data.build_dataset import BuildDataset


class BuildDatasetCellbase(BuildDataset):

    _POST = "post"
    _SCORE = "score"
    _SOURCE = "source"
    _RESULT_FIELD = "result"
    _PHAST_CONS = "phastCons"
    _VARIANT_ID_SEPARATOR = ":"
    _CELLBASE_QUERY_BATCH_SIZE = 200

    _CONSERVATION = "conservation"
    _FUNCTIONAL_SCORES = "functionalScore"
    _GENE_TRAIT_ASSOCIATION = "geneTraitAssociation"
    _VARIANT_TRAIT_ASSOCIATION = "variantTraitAssociation"

    _INCLUDE_LIST = [
        "annotation." + _CONSERVATION,
        "annotation." + _FUNCTIONAL_SCORES,
        "annotation." + _GENE_TRAIT_ASSOCIATION,
        "annotation." + _VARIANT_TRAIT_ASSOCIATION,
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
        self._set_dataset_index_helper_by_attribute("rs_id")
        self._annotate_variation()

    def _annotate_variation(self):
        """
        We build batches of variants to query Cellbase with, so that we don't send too big of a request which
        threatens to send the service down.
        :return:
        """

        list_of_batches = []
        variant_ids_to_query = [
            variant_id for variant_id in self.dataset_index_helper.keys() if variant_id
        ]
        for list_chunk in range(
            0, len(variant_ids_to_query), self._CELLBASE_QUERY_BATCH_SIZE
        ):
            list_of_batches.append(
                variant_ids_to_query[
                    list_chunk: list_chunk + self._CELLBASE_QUERY_BATCH_SIZE
                ]
            )

        # we are putting a hard cap of threads here to not overload Cellbase
        pool = ThreadPool(os.cpu_count())
        pool.map(self._call_cellbase_variation, list_of_batches)

    @renew_access_token
    def _call_cellbase_variation(self, variant_ids_to_query):
        """
        This method takes a list of VariantInfo objects. It will query Cellbase and update the object with the new info.
        Some information of algorithms is only available for GRCh37. As such, we use this method to query Cellbase
        with rs ids, to get the corresponding build 37 variant that have the relevant information.
        :param variants_to_query:
        :return:
        """
        # since we're querying with rs_ids we are able to get relevant information for our build38 variant in a
        # build37 search in Cellbase (since there is more information available at the time)

        # this doesn't yield back, so we have to retrieve all at once. It's only max of 200 cases per query, so
        # it's fine to hold in memory here
        response = self.cellbase_client.search(
            id=variant_ids_to_query, method=self._POST, include=self._INCLUDE_LIST
        )

        # TODO results from cellbase and queries need to have same length, needs fix
        if not response[0][self._RESULT_FIELD] or len(response[0][self._RESULT_FIELD]) != len(variant_ids_to_query):
            raise Exception

        for variant_id_list_index, individual_result in enumerate(
            response[0][self._RESULT_FIELD]
        ):

            # now we fill in all the variant information for these variants
            for variant_info in self.dataset_index_helper[
                variant_ids_to_query[variant_id_list_index]
            ]:
                cellbase_variant_info = individual_result["annotation"]

                # There should always be only one response. However if there are more here, we ignore them
                variant_info.update_object(
                    **{
                        "CADD_scaled_score": self._get_cadd_classification(
                            cellbase_variant_info
                        ),
                        "GERP": self._get_conservation_score_from_source(
                            cellbase_variant_info[self._CONSERVATION], "gerp"
                        ),
                        "phastCons": self._get_conservation_score_from_source(
                            cellbase_variant_info[self._CONSERVATION], "phastCons"
                        ),
                        "phylop": self._get_conservation_score_from_source(
                            cellbase_variant_info[self._CONSERVATION], "phylop"
                        ),
                        "clinVar": self._get_clinvar_classification(
                            cellbase_variant_info
                        ),
                    }
                )

    def _get_conservation_score_from_source(self, conservation, source):
        """
        From a given conservation scores dictionary, this method returns the score for the given source.
        :param conservation:
        :return:
        """
        for conservation_values in conservation:
            if conservation_values[self._SOURCE] == source:
                return conservation_values[self._SCORE]

    def _get_cadd_classification(self, annotation):
        """
        This method extracts the scaled CADD value from a list of functional scores for a given variant.
        :param functional_scores:
        :return:
        """
        if self._FUNCTIONAL_SCORES in annotation:
            for functional_score in annotation[self._FUNCTIONAL_SCORES]:
                if functional_score[self._SOURCE] == "cadd_scaled":
                    return functional_score[self._SCORE]

    def _get_clinvar_classification(self, annotation):
        """
        From the annotation response from Cellbase (if it exists), extract ClinVar information.
        :param annotation:
        :return:
        """
        if self._VARIANT_TRAIT_ASSOCIATION in annotation:
            variant_trait_association = annotation[self._VARIANT_TRAIT_ASSOCIATION]
            if "clinvar" in variant_trait_association:
                # if there are only numbers in the accession,
                # that is the one we want, as opposed to the individual records.
                for record in variant_trait_association["clinvar"]:
                    if record["accession"].isdigit():
                        return record["clinicalSignificance"]
