import re
import numpy as np
import pandas as pd


# FNC: Function which reads _det.m and _dep.m Serpent files
def read_sss_detdep(*filenames):
    """
    filenames: path(s) of file(s) to be read
    """
    output = {}
    for filename in filenames:
        with open(filename, 'r') as file:
            content = file.read()

            # Single variable match
            single_value_matches = re.findall(r'(\w+)\s*=\s*([-+]?\d*\.\d+(?:E[-+]?\d+)?|[-+]?\d+(?:E[-+]?\d+)?);',
                                              content, re.IGNORECASE)

            # Array/matrix match
            array_value_matches = re.findall(r'(\w+)\s*=\s*\[(.*?)\];', content, re.DOTALL)

            # Save single var
            for key, value in single_value_matches:
                value = [int(value)]
                output[key] = np.array([value])

            # Save array/matrix
            for key, values in array_value_matches:
                if key == 'NAMES':
                    row_values = [s.strip() for s in
                                  re.findall(r"'(.*?)'", values.strip())]  # re.findall(r"'(.*?)'", values.strip())
                    output[key] = np.array([row_values]).reshape(-1, )  # , dtype='U')
                else:
                    row_values = re.split(r'%\s*\w+\n|\n', values.strip())  # values.strip().split('\n')

                    all_rows = []
                    for row_str in row_values:
                        if row_str.strip():
                            if key == 'ZAI':
                                row = np.array([int(x) for x in
                                                re.findall(r'[-+]?\d*\.\d+(?:E[-+]?\d+)?|[-+]?\d+(?:E[-+]?\d+)?',
                                                           row_str,
                                                           re.IGNORECASE)])
                            else:
                                row = np.array([float(x) for x in
                                                re.findall(r'[-+]?\d*\.\d+(?:E[-+]?\d+)?|[-+]?\d+(?:E[-+]?\d+)?',
                                                           row_str,
                                                           re.IGNORECASE)])
                            all_rows.append(row)

                    value_matrix = np.array(all_rows)
                    output[key] = value_matrix

    return output


# FNC: Function which reads _res.m Serpent files
def read_sss_res(filename, *search_strings):
    if not search_strings:  # If search_strings is None, gather all unique entries in file
        search_strings = set()
        # unique_strings = set()
        with open(filename, 'r') as file:
            for line in file:
                identifier = re.match(r'(\w+)(\s+\(idx,)', line)  # Find entry at line beginning
                if identifier:
                    # unique_strings.add(identifier.group(1))
                    search_strings.add(identifier.group(1))
        search_strings = sorted(search_strings)
        # search_strings = unique_strings

    result = {string: [] for string in search_strings}
    with open(filename, 'r') as file:
        lines = file.readlines()

    idx = -1  # 0

    # Search str in file, get values after "=", append to results
    for line in lines:
        if line.startswith('% Increase counter'):
            idx += 1
        else:
            # print(line)
            for search_string in search_strings:
                if line.startswith(search_string):
                    if re.match(r'(\w+)', line).group() == search_string:  #line.startswith(search_string):  #
                        values_line = line.split('=')[1].strip().rstrip(';')
                        # Check if float, string or list
                        if values_line.startswith('['):
                            values = re.findall(r'[\d\.\-\+E]+', values_line)
                            if len(result[search_string]) == idx:
                                result[search_string].append(list(map(float, values)))
                            else:
                                result[search_string].append(list(np.zeros_like(list(map(float, values)))))
                                result[search_string].append(list(map(float, values)))
                                # result[search_string][idx-1].extend(list(map(float, values)))
                        else:
                            try:
                                value = [float(values_line)]
                                # if len(result[search_string]) == idx:
                                result[search_string].append(value)
                                # else:
                                #     result[search_string][idx-1].extend(value)
                            except ValueError:
                                value = values_line[1:-2]
                                # if len(result[search_string]) < idx:
                                result[search_string].append(value)
                                # else:
                                #     result[search_string][idx-1].extend(value)

    for search_string in search_strings:
        result[search_string] = np.array(result[search_string])

    return result


# FNC: Function which reads .out Serpent files
def read_out_file(fn):
    """
    fn: path of file to be read
    """

    with open(fn, 'r') as ff:  # +'.out'
        lines_out = ff.readlines()

    dict_mat = {}
    for ii, line in enumerate(lines_out):
        if 'Material "' in line:
            # print(ii, line)
            mat = line.split('"')[1]
            try:
                dict_mat[mat] = {
                    'adens': float(lines_out[ii + 5].split(' ')[-2]),
                    'mdens': float(lines_out[ii + 6].split(' ')[-2]),
                    'vol': float(lines_out[ii + 7].split(' ')[-2]),
                    'n_nuc': int(lines_out[ii + 11].split(' ')[2])
                }
            except ValueError:
                print('row ', ii)
            df_tmp = pd.read_table(fn,
                                   sep='\s+',
                                   skiprows=ii + 20,
                                   # nrows=counter,
                                   nrows=dict_mat[mat]['n_nuc'],
                                   names=['Nuclide', 'aweight', 'temp', 'adens', 'afrac', 'mfrac'],
                                   #  index_col=False,
                                   index_col=0,
                                   usecols=[0, 1, 2, 3, 4, 5],
                                   #  converters={'mfrac': lambda x: re.sub(r'\D+', '', x)}
                                   )
            df_tmp.index = [nn.split('.')[0] for nn in list(df_tmp.index)]
            dict_mat[mat]['comp'] = df_tmp

    return lines_out, dict_mat
