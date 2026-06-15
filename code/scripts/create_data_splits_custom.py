import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import shuffle
import math
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

FILE = 'swissprot_full.tsv'
SEED = 10297362
df = pd.read_csv(FILE, sep='\t', header=0)

n_splits=5
folds = []

orphans_used = []
for i in range(n_splits):
    df = shuffle(df, random_state=SEED + i)  # seed per fold for reproducibility (Bug F)
    extra_df = pd.DataFrame(columns=df.columns)
    train = pd.DataFrame(columns=df.columns)
    valid = pd.DataFrame(columns=df.columns)
    test = pd.DataFrame(columns=df.columns)
    orphans = pd.DataFrame(columns=df.columns)

    #parse out test2 "out of distribution"
    for ec_number in df['EC number'].unique():
        ec_df = df[df['EC number'] == ec_number]
        grouped = ec_df.groupby('UniRef50')
        num_groups = grouped.ngroups
    
        counter = 'Test'
        if num_groups < 3:
            extra_df = pd.concat([extra_df, ec_df])
            continue

        if num_groups >= 10:
            num_test = math.ceil(0.1 * num_groups)
            num_valid = math.ceil(0.1 * num_groups)
            num_train = num_groups - num_test - num_valid
        elif num_groups >= 6 and num_groups <10 :
            num_test = 2
            num_valid = 2
            num_train = num_groups - num_test - num_valid
        else:
            num_test = 1
            num_valid = 1
            num_train = num_groups - num_test - num_valid

        for group_label, group_df in grouped:
            if num_test > 0 and counter == 'Test':
                test = pd.concat([test, group_df])
                num_test -= 1
                if num_test == 0:
                    counter = 'Valid'
            elif num_valid > 0 and counter == 'Valid':
                valid = pd.concat([valid, group_df])
                num_valid -= 1
                if num_valid == 0:
                    counter = 'Train'
            else:
                train = pd.concat([train, group_df])

    # Group the dataframe by 'EC'
    df2 = extra_df.groupby('EC number')
    # Count the occurrences of each 'EC' group
    counts_df = df2.size().reset_index(name='count')
    
    # Iterate over each group
    for ec_number, group in df2:
        count = counts_df[counts_df['EC number'] == ec_number]['count'].values[0]
        if count >= 3:
            test = pd.concat([test, group.iloc[0:1]])
            valid = pd.concat([valid, group.iloc[1:2]])
            train = pd.concat([train, group.iloc[2:]])
        elif count == 2:
            train = pd.concat([train, group.iloc[:1]])
            test = pd.concat([test, group.iloc[1:]])
        else:
            orphans = pd.concat([orphans, group])

    #split 20% orphans into test2 and 80% into training
    n_sample = int(len(orphans)/n_splits)
    to_sample = orphans[~orphans.index.isin(orphans_used)]

    if n_sample > len(to_sample):
        n_sample = len(to_sample)

    test2 = to_sample.sample(n=n_sample, random_state=SEED + i)  # Bug F: seed the sample
    train = pd.concat([train, orphans[~orphans.index.isin(test2.index)]])

    orphans_used.extend(list(test2.index)) #dont reuse orphans from previous splits

    print('Split: ', i+1)
    print('Train:', train.shape)
    print('Valid:', valid.shape)
    print('Test1:', test.shape)
    print('Test2:', test2.shape)
    print()

    train_ids = train.index.to_list()
    valid_ids = valid.index.to_list()
    test_ids = test.index.to_list()
    test2_ids = test2.index.to_list()

    folds.append((train_ids, valid_ids, test_ids, test2_ids)) #for plotting

    train.to_csv('./custom/split_{}/train.csv'.format(i+1), sep='\t', index=False)
    valid.to_csv('./custom/split_{}/valid.csv'.format(i+1), sep='\t', index=False)
    test.to_csv('./custom/split_{}/test1.csv'.format(i+1), sep='\t', index=False)
    test2.to_csv('./custom/split_{}/test2.csv'.format(i+1), sep='\t', index=False)  # Bug A: was `orphans` (full set), leaking ~80% of test-2 into train

#final dataset will be stratify k-fold
le = LabelEncoder().fit(df['EC number'])
#pool_df = df[~df.index.isin(orphans.index.to_list())]

X = df['Entry']
y_EC = df['EC number']
y = le.transform(y_EC)
n_classes = len(le.classes_)

cmap_data = plt.cm.Paired
cmap_cv = plt.cm.coolwarm
#https://scikit-learn.org/stable/auto_examples/model_selection/plot_cv_indices.html#sphx-glr-auto-examples-model-selection-plot-cv-indices-py
def plot_cv_indices(folds, X, y, group, ax, n_splits=5, lw=10, use_groups=False):
    
    """Create a sample plot for indices of a cross-validation object."""
    groups = group if use_groups else None
    # Generate the training/testing visualizations for each CV split
    for ii, (train, valid, test, test2) in enumerate(folds):
        # Fill in indices with the training/test groups
        indices = np.array([np.nan] * len(X))
        indices[train] = 0.0
        indices[valid] = 0.4
        indices[test] = 0.8
        indices[test2] = 1.0
    
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
    yticklabels = list(range(1, n_splits+1)) + ["EC Number", "UniRef50 Cluster"]
    ax.set(
        yticks=np.arange(n_splits + 2) + 0.5,
        yticklabels=yticklabels,
        xlabel="Sample index",
        ylabel="CV iteration",
        ylim=[n_splits + 2.2, -0.2],
        xlim=[0, len(X)],
    )
    ax.set_title("{}".format('EC-Stratified UniRef50-KFold'), fontsize=15)
    return ax

fig, ax = plt.subplots(figsize=(6, 3))
plot_cv_indices(folds, X, y, LabelEncoder().fit_transform(df['UniRef50']), ax)
ax.legend(
    [Patch(color=cmap_cv(0.0)), Patch(color=cmap_cv(0.4)), Patch(color=cmap_cv(0.8)), Patch(color=cmap_cv(1.0))],
    ["Train", "Validation", 'Test 1', 'Test 2'],
    loc=(1.02, 0.5),
)
# Make the legend fit
plt.tight_layout()
fig.subplots_adjust(right=0.7)
plt.savefig('./custom/train_splits.png')
