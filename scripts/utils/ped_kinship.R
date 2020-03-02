library('GENLIB')

ped_data <- read.table(
        '~/project/pedigree_msp/data/Luke/Genizon_4149gen_nov2019.txt',
        header=1)

head(ped_data)
ped_data['ind'] <- ped_data['id']
pedigree <- gen.genealogy(ped_data)

probands <- gen.pro(pedigree)

kinship <- gen.phi(pedigree, probands[1:10])
kinship
