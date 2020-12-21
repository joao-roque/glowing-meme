class VariantEntryInfo:

    VARIANT_INFO_VALUES = [
        "id",
        "chromosome",
        "start",
        "end",
        "alt",
        "ref",
        "assembly",
        "case_id",
        "rs_id",
        "age",
        "sex",
        "zygosity_proband",
        "zygosity_mother",
        "zygosity_father",
        "tier",
        "mode_of_inheritance",
        "consequence_type",
        "biotypes",
        "population_frequency",
        "CADD_score",
        "type",
        "clinVar",
        "DisGeNET",
        "PhastCons",
        "phylop",
        "GERP",
        "program",
        "mother_ethnic_origin",
        "father_ethnic_origin",
        "segregation_pattern",
        "penetrance",
        "gel_variant_acmg_classification",
        "case_solved_family",
        "phenotypes_solved",
        "actionability",
        "interpretation_message",
        "dict_extra_scores",
        "reported_outcome",
    ]

    def __init__(self, **kwargs):
        """
        This object holds the variant entry info and can be set with a dictionary
        :param kwargs:
        """
        # setting variables as defined in _VARIANT_INFO_VALUES
        # this allows flexibility to add more values in the future
        for key in self.VARIANT_INFO_VALUES:
            if key in kwargs:
                setattr(self, key, kwargs[key])
            else:
                setattr(self, key, None)

    def update_object(self, **kwargs):
        """
        This method updates attributes given in dict kwargs.
        :return:
        """
        for key in self.VARIANT_INFO_VALUES:
            if key in kwargs:
                setattr(self, key, kwargs[key])

    def __iter__(self):
        """
        This method returns an ordered iterator of this class's attributes as per VARIANT_INFO_VALUES.
        :return:
        """
        return iter([value for attr, value in self.__dict__.items() if attr != 'VARIANT_INFO_VALUES'])

