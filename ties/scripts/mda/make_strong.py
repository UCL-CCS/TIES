import parmed as pmd

t = pmd.load_file('sys_solv.top')

# find the angle
for ang in t.angles:
	conf = ['CL8', 'BR1']
	ang_conf = [ang.atom1.name, ang.atom3.name]
	if ang_conf == conf or ang_conf == conf[::-1]:
		print('found ang', ang)
		chosen = ang
	pass