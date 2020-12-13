import time
from glowingmeme.build_data.build_dataset import BuildDatasetCVA, BuildDatasetCipapi


start_time = time.time()
bd_cva = BuildDatasetCVA()
main_dataset = bd_cva.build_dataset()

bd_cipapi = BuildDatasetCipapi(main_dataset)
main_dataset = bd_cipapi.build_dataset()
print("--- %s seconds ---" % (time.time() - start_time))

x = 1
