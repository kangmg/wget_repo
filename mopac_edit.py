with open("mopac.out", "r+") as f:
    lines = f.readlines()
    f.seek(0)
    for line in lines:
        if "= FINAL HEAT OF FORMATION" not in line:
            f.write(line)
    f.truncate()
