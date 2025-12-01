# params.py

fractions = ["Adipocytes", "Immune", "Stromal", "Epithelial"]
subtype_order = ["Immunoreactive", "Differentiated", "Mesenchymal", "Proliferative"]
subtype_map = {
    "IMR_consensus": "Immunoreactive",
    "DIF_consensus": "Differentiated",
    "MES_consensus": "Mesenchymal",
    "PRO_consensus": "Proliferative",
}
cont_cols = ["Age", "BMI"]
custom_colors = ["green", "#FFD700", "darkorange", "red"]
colors = ["peachpuff", "orange", "tomato", "salmon"]
levels = subtype_order
no_adj_samples = [130131, 130138, 161063, 161090, 170110, 190001, 190010]
