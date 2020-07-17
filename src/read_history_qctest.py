def decode_qctest(hex_code):
    # hex to numeric
    num = int(hex_code, 16)
    # list to save test number in
    tests = []
    for i in range(26,0,-1):
        qc_binary_id = 2**i
        if qc_binary_id <= num:
            num -= qc_binary_id
            tests.append(i)
    
    return tests[::-1], num

def print_result(tests, remainder):
    print('Tests performed: ', end='')
    for t in tests[:-1]:
        print('{:d}, '.format(t), end='')
    print('{:d}'.format(tests[-1]))
    print('Remainder (should be 0): {:d}'.format(remainder))

# get test numbers and print them out
print('FFBFE:')
print_result(*decode_qctest('FFBFE'))
print('\n4000:')
print_result(*decode_qctest('4000'))