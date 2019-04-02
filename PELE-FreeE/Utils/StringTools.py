import re


def natural_sort(l):
    def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
        return [int(text) if text.isdigit() else text.lower()
                for text in _nsre.split(s)]

    return sorted(l, key=natural_sort_key)
