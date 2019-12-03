
for name in ['mini', 'small', 'medium', 'large', 'extra']:
    time_test = []
    for i in range(1, 22):
        time_test.append([])
        for k in range(1, 21):
            with open('./{name}_out/{name}{i}_{k}.txt'.format(name=name, i=i, k=k), 'r') as f:
                for line in f:
                    if 'Time in seconds = ' in line:
                        x = eval(line.split(' = ')[1])
                        time_test[i - 1].append(x)
                        break
        time_test[i - 1].append(sum(time_test[i - 1]) / 20)

    for i in range(0, 21):
        with open('results/{name}{i}.txt'.format(name=name, i=i), 'w') as f:
            res = ""
            for j in range(0, 21):
                res += str(8 * (j + 1)) + '\t' + str(time_test[j][i]) + '\n'
            f.write(res)