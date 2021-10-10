import pandas as pd
from sklearn.decomposition import PCA

#set border length 
border_length = 4

# generate list of shapes
all_shapes = ['Stagger', 'Rise', 'Opening', 'Buckle', 'MGW', 'Tilt', 'HelT', 'Roll', 'Shear', 'Slide', 'Stretch', 'ProT', 'Shift']

# generate whole set
shape_fimo = pd.read_csv("shape_fimo_forward.csv")

core_motif_len = len(shape_fimo["matched_sequence"].iloc[0])

header = []
for shape in all_shapes:
  for i in range((30 - border_length),(30 + core_motif_len + border_length)):
    header.append(shape + "_" + str(i))

X_shapes = shape_fimo[header]


#PCA
pca = PCA(n_components=10)
principal_components = pca.fit_transform(X_shapes.to_numpy())
principal_df = pd.DataFrame(data = principal_components, columns = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'])

principal_df.to_csv("dimensionally_reduced_features.tsv", sep="\t", index=False)
