from glowingmeme.build_data.variant_entry_info import VariantEntryInfo


variant_info_list = [VariantEntryInfo(**{"program": "100k"})]

dataset_by_key = {}
dataset_key = "program"
if dataset_key in VariantEntryInfo.VARIANT_INFO_VALUES:
    for variant_info_object in variant_info_list:
        if vars(variant_info_object)[dataset_key] in dataset_by_key:
            dataset_by_key[dataset_key].append(variant_info_object)
        else:
            dataset_by_key[dataset_key] = [variant_info_object]

print(dataset_by_key["program"][0].age)
dataset_by_key["program"][0].age = 1
print(variant_info_list[0].age)
