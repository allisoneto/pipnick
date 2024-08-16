# Conditions used to sort data by a category 
# Uses range of image numbers (d1078.fits = 78)

conditions_06_26 = [(1.375, (1065, 1074)),
                    (1.624, (1022, 1031)),
                    (1.625, (1088, 10105)),
                    (1.875, (1033, 1042)),
                    (2.625, (1043, 1053)),
                    (3.375, (1054, 1064)),
                    ]

conditions_06_24 = [(1.375, (1053, 1060)),
                    (1.625, (1001, 1052)),
                    (1.625, (1088, 1105))
                    ]

conditions = {'06-26': conditions_06_26, '06-24': conditions_06_24}