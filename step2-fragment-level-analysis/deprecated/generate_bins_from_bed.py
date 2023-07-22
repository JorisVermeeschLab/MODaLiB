regions = list()
with open('atlas.bed', 'r') as f:
    for line in f.readlines():
        elements = line.rstrip().split()
        chr_name = elements[0]
        if chr_name[:3] == 'chr':
            chr_name = chr_name[3:]
        chr_name = int(chr_name)
        start = int(elements[1])
        end = int(elements[2])
        regions.append((chr_name, start, end))


regions = list(sorted(regions, key=lambda x: x[0] * 1e+9 + x[1]))


with open('preselected-bins.txt', 'w') as f:
    for region in regions:
        f.write('%i, %i, %i\n' % tuple(region))



with open('preselected-bins-small.txt', 'w') as f:
    for region in regions:
        chr_name, start, end = region
        f.write('%i, %i, %i\n' % (chr_name, start + 200, end - 200))
