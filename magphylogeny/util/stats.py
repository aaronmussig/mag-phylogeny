def weighted_jaccard(vec_a, vec_b):
    minima, maxima = list(), list()
    for a, b in zip(vec_a, vec_b):
        minima.append(min(a, b))
        maxima.append(max(a, b))
    try:
        return sum(minima) / sum(maxima)
    except ZeroDivisionError:
        return 0.0
