def read_calc_freq(file_name):
    with open(file_name) as input_file:
        lines = []
        # input_file.readline(2)
        for i, line in enumerate(input_file.readlines()[2:]):
            if line != '\n' and i > 1:
                lines.append(line)
            elif i > 1:
                break
        result = [float(freq) for line in lines for freq in line.strip('\n').split()]
    return sorted(result)


def read_exp_freq(file_name):
    with open(file_name) as input_file:
        result = []
        for line in input_file.readlines():
            result.append(float(line.strip('\n').split()[0]))
    return sorted(result)


def compare_freq(structure_name, path='../../structures/', print_log=True):
    path = path + structure_name + '/'
    exp_freq_path = path + 'intencities'
    calc_freq_path = path + 'skra.ved'

    exp_freq = read_exp_freq(exp_freq_path)
    calc_freq = read_calc_freq(calc_freq_path)

    max_delta = abs(calc_freq[0] - exp_freq[0])
    for i in range(len(calc_freq)):
        if print_log:
            print(i, round(abs(calc_freq[i] - exp_freq[i]), 2))
        max_delta = max(max_delta, abs(calc_freq[i] - exp_freq[i]))
    if print_log:
        print('max delta', round(max_delta, 2))
    return max_delta


if __name__ == '__main__':
    compare_freq('gketoohw12')