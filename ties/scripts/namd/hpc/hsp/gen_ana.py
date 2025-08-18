# tyk2
trans = ["l15_l10", "l15_l16", "l15_l6", "l1_l10", "l1_l3"]
# mcl?
# l12_l35  l16_l34  l17_l9  l18_l39  l32_l38  l38_l42  l39_l42  l6_l41
# thrombin
trans = ["l1_l8", "l1_l9", "l3_l5", "l4_l10", "l5_l6"]
# ptpb1
trans = ["l13_l20", "l3_l23", "l3_l7", "l4_l22", "l8_l14"]
# cdk2
trans = ["l1q_l21", "l1q_l26", "l20_l21"]

counter = 1
for tran in trans:
    s = open("run_script_template.sh").read()
    open(f"run_script_{counter}.sh", "w").write(s.format(transformation=tran))
    counter += 1
