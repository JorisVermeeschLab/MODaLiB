def pair(b1, b2):
    n = len(b1)
    m = len(b2)

    s1 = ''
    s2 = ''

    i = 0
    j = 0

    while (i < n) or (j < m):

        if i < n:
            min1 = b1[i]
        else:
            min1 = float('inf')
        if j < m:
            min2 = b2[j]
        else:
            min2 = float('inf')
        min_both = min(min1, min2)

        use1 = (min_both == min1)
        use2 = (min_both == min2)


        if use1:
            s1 += str(b1[i]) + ' '
        else:
            s1 += '-- '
        if use2:
            s2 += str(b2[j]) + ' '
        else:
            s2 += '-- '


        if use1:
            i += 1
        if use2:
            j += 1

    print(s1)
    print(s2)



bases1 = [34, 35, 36, 36, 36, 36, 36, 37, 38, 46, 47, 48]
bases2 = [47, 48, 49, 50, 51, 52, 53, 54]


pair(bases1, bases2)
