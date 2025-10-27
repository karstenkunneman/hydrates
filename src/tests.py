import model

#1-Component T-Defined Fresh Water
assert round(model.KlaudaSandler2003([1], [1], "T", 280, None).pressure/1E5, 1) == 51.5 #bar

#5-Component T-Defined Fresh Water
assert round(model.KlaudaSandler2003([1,2,3,7,9], [0.932, 0.0425, 0.0161, 0.0051, 0.0043], "T", 277.7, None).pressure/1E5, 1) == 16.3 #bar

#1-Component P-Defined Fresh Water
assert round(model.KlaudaSandler2003([3], [1], "P", None, 2E5).temperature, 1) == 274.1 #K

#3-Component P-Defined Fresh Water
assert round(model.KlaudaSandler2003([1,2,3], [0.9196,0.0513,0.0291], "P", None, 10.7E5).temperature, 1) == 275.8 #K