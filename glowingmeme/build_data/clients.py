import os
import yaml
from pyark.cva_client import CvaClient
from pycellbase.cbclient import CellBaseClient
from pycellbase.cbconfig import ConfigClient
from pycipapi.cipapi_client import CipApiClient


class Clients:
    """
    Connect to production clients using the GEL_CREDENTIALS as an env variable
    """

    PROD_CVA_HOST = "https://bio-prod-cva.gel.zone"
    PROD_CIPAPI_HOST = "https://cipapi-prod.gel.zone"
    PROD_CELLBASE_HOST = "https://cellbase.gel.zone/"

    def __init__(self, cipapi_host=None, cva_host=None, cellbase_host=None):
        self._credentials = self._get_credentials()
        self.cva_host = cva_host if cva_host else self.PROD_CVA_HOST
        self.cipapi_host = cipapi_host if cipapi_host else self.PROD_CIPAPI_HOST
        self.cellbase_host = cellbase_host if cellbase_host else self.PROD_CELLBASE_HOST

    def _get_credentials(self):
        """
        Build credentials from GEL_CREDENTIALS env file.
        :return:
        """
        credentials = {entry['name']: entry for entry in yaml.load(open(os.getenv('GEL_CREDENTIALS')),
                                                                   Loader=yaml.FullLoader)}
        return credentials

    def get_cipapi_client(self):
        """
        Get and login to cipapi client
        :return:
        """
        return CipApiClient(url_base=self.cipapi_host, user=self._credentials['cip_api_prod']['username'],
                            password=self._credentials['cip_api_prod']['password'], retries=8)

    def get_cellbase_client(self):
        """
        Get and login to cellbase variant client
        :return:
        """
        cellbase_configuration = {"species": "hsapiens",
                                  "version": "v4",
                                  "rest": {"hosts": [self.PROD_CELLBASE_HOST]}}
        cellbase_client = CellBaseClient(ConfigClient(cellbase_configuration))
        cellbase_variant_client = cellbase_client.get_variant_client()

        return cellbase_variant_client

    def get_cva_client(self):
        """
        Get and login to CVA client
        :return:
        """
        return CvaClient(url_base=self.PROD_CVA_HOST, user=self._credentials['cva_prod']['username'],
                         password=self._credentials['cva_prod']['password'])

    def get_all_clients(self):
        """
        Get cipapi, cellbase and cva production clients.
        :return: cipapi_client, cellbase_client, cva_client
        """
        return self.get_cipapi_client(), self.get_cellbase_client(),  self.get_cva_client()
