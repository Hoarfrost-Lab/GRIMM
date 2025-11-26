import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold, KFold, train_test_split
from sklearn.preprocessing import LabelEncoder
import math
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

FILE = 'swissprot_full.tsv'
SEED = 10297362
df = pd.read_csv(FILE, sep='\t', header=0)

extra_df = pd.DataFrame(columns=df.columns)
pool_df = pd.DataFrame(columns=df.columns)
orphans = pd.DataFrame(columns=df.columns)

#parse out test2 "out of distribution"
for ec_number in df['EC number'].unique():
    ec_df = df[df['EC number'] == ec_number]
    grouped = ec_df.groupby('UniRef50')
    num_groups = grouped.ngroups
    
    if num_groups < 3:
        extra_df = pd.concat([extra_df, ec_df])
        continue

    for group_label, group_df in grouped:
        pool_df = pd.concat([pool_df, group_df])

# Group the dataframe by 'EC'
df2 = extra_df.groupby('EC number')
# Count the occurrences of each 'EC' group
counts_df = df2.size().reset_index(name='count')
    
# Iterate over each group
for ec_number, group in df2:
    count = counts_df[counts_df['EC number'] == ec_number]['count'].values[0]
    if count >= 2:
        pool_df = pd.concat([pool_df, group])
    else:
        orphans = pd.concat([orphans, group])

print('Pool:', pool_df.shape)
print('Orphans:', orphans.shape)

#final dataset will be stratify k-fold
le = LabelEncoder().fit(df['EC number'])

X = pool_df['Entry']
y_EC = pool_df['EC number']
y = le.transform(y_EC)
n_classes = len(le.classes_)
group = pool_df['UniRef50']


skf = StratifiedKFold(n_splits=5)
folds = skf.split(X, y)

kf = KFold(n_splits=5)
orphan_folds = kf.split(orphans)

for i, ((train_index, test_index), (o_train_index, o_test_index)) in enumerate(zip(folds, orphan_folds)):
    X_train = pool_df.iloc[list(train_index)] #80% data reserved for training

    counts = Counter(y[list(test_index)])
    modified_y = [_y if counts[_y] > 1 else n_classes for _y in y[list(test_index)]] #replace validation singltons with a new "class" so they get evenly split

    X_valid, X_test1 = train_test_split(pool_df.iloc[list(test_index)], test_size=0.5, random_state=SEED, stratify=modified_y)
    
    X_train = pd.concat([X_train, orphans.iloc[list(o_train_index)]]) #adding some orphans to the training data
    X_test2 = orphans.iloc[list(o_test_index)] #20% go into test2 in each fold

    print('Split:', i+1)
    print('Train:', X_train.shape)
    print('Valid:', X_valid.shape)
    print('Test1:', X_test1.shape)
    print('Test2:', X_test2.shape)
    print()

    X_train.to_csv('./scikit/split_{}/train.csv'.format(i+1), index=False, sep='\t')
    X_valid.to_csv('./scikit/split_{}/valid.csv'.format(i+1), index=False, sep='\t')
    X_test1.to_csv('./scikit/split_{}/test1.csv'.format(i+1), index=False, sep='\t')
    X_test2.to_csv('./scikit/split_{}/test2.csv'.format(i+1), index=False, sep='\t')

cmap_data = plt.cm.Paired
cmap_cv = plt.cm.coolwarm
#https://scikit-learn.org/stable/auto_examples/model_selection/plot_cv_indices.html#sphx-glr-auto-examples-model-selection-plot-cv-indices-py
def plot_cv_indices(cv, X, y, group, ax, n_splits=5, lw=10):
    
    """Create a sample plot for indices of a cross-validation object."""
    use_groups = "Group" in type(cv).__name__
    groups = group if use_groups else None
    # Generate the training/testing visualizations for each CV split
    for ii, (tr, tt) in enumerate(cv.split(X, y)):
        # Fill in indices with the training/test groups
        indices = np.array([np.nan] * len(X))
        indices[tt] = 1
        indices[tr] = 0

        # Visualize the results
        ax.scatter(
            range(len(indices)),
            [ii + 0.5] * len(indices),
            c=indices,
            marker="_",
            lw=lw,
            cmap=cmap_cv,
            vmin=-0.2,
            vmax=1.2,
        )

    # Plot the data classes and groups at the end
    ax.scatter(
        range(len(X)), [ii + 1.5] * len(X), c=y, marker="_", lw=lw, cmap=cmap_data
    )

    ax.scatter(
        range(len(X)), [ii + 2.5] * len(X), c=group, marker="_", lw=lw, cmap=cmap_data
    )

    # Formatting
    yticklabels = list(range(n_splits)) + ["class", "group"]
    ax.set(
        yticks=np.arange(n_splits + 2) + 0.5,
        yticklabels=yticklabels,
        xlabel="Sample index",
        ylabel="CV iteration",
        ylim=[n_splits + 2.2, -0.2],
        xlim=[0, 100],
    )
    ax.set_title("{}".format('SciKit StatifiedKFold'), fontsize=15)
    return ax

fig, ax = plt.subplots(figsize=(6, 3))
plot_cv_indices(skf, X, y, LabelEncoder().fit_transform(pool_df['UniRef50']), ax)
ax.legend(
    [Patch(color=cmap_cv(0.8)), Patch(color=cmap_cv(0.02))],
    ["Testing set", "Training set"],
    loc=(1.02, 0.8),
)
# Make the legend fit
plt.tight_layout()
fig.subplots_adjust(right=0.7)
plt.savefig('./scikit/train_splits.png')
