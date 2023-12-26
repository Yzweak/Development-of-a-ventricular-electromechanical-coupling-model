import functools


def compare(a, b):
    if a != b:
        return -1 if a > b else 1

num_list = [4, 2, 8, -9, 1, -3]
sorted_num_list = sorted(num_list,key=functools.cmp_to_key(compare))
print(sorted_num_list)
