def convert_h_to_ru(file_name='input.dd2'):
    # h_c_compositions = ['HC', 'CH', 'CCOH']
    # replacement = ['RuC', 'CRu', 'CCORu']
    ru_numbers = [x for x in range(54, 72)]
    with open(file_name) as input_file:
        with open('output.dd2', 'w') as output_file:
            for line in input_file.readlines():
                line_list = line.split()
                if 12 > len(line_list) > 3:
                    name_list = [x for x in line_list[-2]]
                    for x in ru_numbers:
                        if (line_list[0] in ['s', 'v', 'k']
                            and str(x) in line_list[4:]) \
                                or ((line_list[0] == '1.00' or line_list[0] == '-1.00') and str(x) in line_list[1:]):
                            x_index = line_list.index(str(x))
                            index_in_name = len(name_list) + (x_index - len(line_list) + 2)
                            name_list[index_in_name] = 'Ru'
                    line = line.replace(line_list[-2], ''.join(name_list))
                output_file.write(line)


if __name__ == '__main__':
    convert_h_to_ru()

