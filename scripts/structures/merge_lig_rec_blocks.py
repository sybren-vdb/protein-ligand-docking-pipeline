import sys
import gemmi
import os

def chain_names(model): return {ch.name for ch in model}

for fn in sys.argv[1:]:
	doc = gemmi.cif.read_file(fn)
	rec = gemmi.make_structure_from_block(doc[0])
	model = rec[0]

	if len(doc) > 1:
		lig_block = doc[1]
		if lig_block.find_loop('__atom_site.') is not None:
			lig = gemmi.make_structure_from_block(lig_block)
			used = chain_names(model)
			for ch in lig[0]:
				new = ch.clone()
				base = new.name or 'X'
				name = base
				i = 1
				while name in used:
					name = f'{base}{i}'
					i += 1
				new.name = name
				model.add_chain(new)
	out = rec.make_mmcif_document()
	out.write_file(os.path.splitext(fn)[0] + "_merged.cif")
