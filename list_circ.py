circdir='/atium/Data/Nanopore/Analysis/170104_reassemble/circlate/'
outdir='/home/yfan/Dropbox/Lab/oxa/'

circinfo=[]
for samp in [1,2,4,6,7,8,9,10,12] :
    piloncirc=circdir+str(samp)+'.pilon_circ/04.merge.circularise.log'
    with open(piloncirc) as f:
        content=[i.strip('\n').split('\t') for i in f]
    for i in content[1:]:
        i.insert(1, str(samp))
        circinfo.append(i[1:])

circinfo.insert(0, ['samp', 'contig', 'repeat', 'circ_by_nucmer', 'circ_by_spades', 'circ'])

with open(outdir+'circ_summary.csv', 'w') as f:
    for i in circinfo:
        f.write(','.join(i)+'\n')
    
