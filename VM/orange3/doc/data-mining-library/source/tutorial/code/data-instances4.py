import Orange
from collections.abc import Counter

data = Orange.data.Table("lenses")
print(Counter(str(d.get_class()) for d in data))
