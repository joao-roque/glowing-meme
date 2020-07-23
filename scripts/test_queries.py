from glowingmeme.build_data.clients import Clients


cipapi, cellbase, cva = Clients().get_all_clients()

evidence = cva.evidences()
evidences = evidence.get_evidences(source="CADD")
x = 1
