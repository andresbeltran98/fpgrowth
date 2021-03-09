from fptree import *

dataset = [['bread', 'milk'],
           ['bread', 'diaper', 'beer', 'eggs'],
           ['milk', 'diaper', 'beer', 'coke'],
           ['bread', 'milk', 'diaper', 'beer'],
           ['bread', 'milk', 'diaper', 'coke']]

SUP = 0.4
MIN_CONF = 0.6
SUP_COUNT = round(SUP * len(dataset))

data = data_to_trans(dataset)

# Create header table
ordered = header_table(data, SUP_COUNT)

# Sort transactions based on header table
sort_transactions(data, ordered)

# Construct FPTree
t = FPTree()
t.insert(data, ordered)

# Perform fp_growth algorithm
freq_itemset = []
fp_growth(t, ordered, [], freq_itemset, SUP_COUNT)
print('FPTree Algorithm')
print('Frequent items: ', len(freq_itemset))
for x in freq_itemset:
    print(x)

# Create association rules
print('\nAssociation Rules:')
r = association_rules(freq_itemset, dataset, MIN_CONF)
for rule in r:
    print(rule)