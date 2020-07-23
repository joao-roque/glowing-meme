

class BuildDataset:

    def __init__(self, cipapi, cva):
        """
        Clients used to query services.
        :param cipapi:
        :param cva:
        """
        self.cipapi = cipapi
        self.cva = cva

    def get_variants_used_in_diagnosis(self):
        """
        This methods queries Genomics England services
        to gather variants that were used in the diagnosis by the clinician.
        :return:
        """
        pass
