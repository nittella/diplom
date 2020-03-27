def get_atom_type(struct_name):
    path = '../structures/' + struct_name + '/coords'
    with open(path, 'r') as input_file:
        result = '"' + struct_name + '": {\n'
        lines = [line.strip('\n').split() for line in input_file.readlines()]
        for n, line in enumerate(lines):
            if line[0] == '1':
                result += '\t' + str(n + 1) + ': ("water", "H"),\n'
            elif line[0] == '6':
                result += '\t' + str(n + 1) + ': ("carbon", "C"),\n'
            elif line[0] == '8':
                result += '\t' + str(n + 1) + ': ("water", "O"),\n'
            else:
                print('Uknown atom', line[0])
        result += '}'
        result = result.replace('"', "'")

        return result


if __name__ == '__main__':
    print(get_atom_type('gketoohw1'))
