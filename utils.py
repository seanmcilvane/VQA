def prob_dist(counts, shots):
    output_distr = [v / shots for v in counts.values()]
    if len(output_distr) == 1:
        output_distr.append(1 - output_distr[0])
    return output_distr