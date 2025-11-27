#!/usr/bin/env python3

import click
import genomerator as gen

def zero_data (feature):
	feature.data['count'] = 0
	return feature

def increment_count (region,feature):
	try:
		this_count = int(feature.data[3])
	except ValueError: # if there is no count (score) specified for this feature, treat it as a single count
		this_count = 1
	region.data.__setitem__('count', region.data['count'] + this_count) # increment by count

@click.command()
@click.option('-p', '--position', type = int, default = 0, help = 'position offset of counted base within motif')
@click.option('-z', '--zeroes', is_flag = True, help = 'output regions with zero count')
@click.option('-q', '--quiet', is_flag = True, help = 'don\'t show progress bar')
@click.argument('reference_file', type = click.File('r'))
@click.argument('region_bed', type = click.File('r'))
@click.argument('count_bed', required = False, type = click.File('r'), default = click.open_file('-', 'r'))
@click.argument('out_tsv', required = False, type = click.File('w'), default = click.open_file('-', 'w'))
def region_counts (region_bed, count_bed, out_tsv, reference_file, position, zeroes, quiet):
	'''
	Given a BED file of genome regions, a BED file of hits per feature (e.g. restriction motif site), and a list of reference lengths, count the number of total hits in each genome region.
	Each feature is counted at the single specified position relative to its start, in the same orientation:
	e.g. a 6-base feature from chr1:1001 to chr1:1006 with offset "-p 1" would be counted at position 1002,
	but if it is reverse oriented from chr1:1006 to chr1:1001 it would be counted at position 1005.
	Result is written as two-column TSV with region name and hit count.
	'''
	
	references = gen.read_reference_lengths(reference_file)
	reference_names = list(references.keys())
	reference_lengths = list(references.values())
	region_stream = gen.BedStream(region_bed, assert_sorted = True, references = references.keys())
	count_stream = gen.BedStream(count_bed, assert_sorted = True, parse = False, references = references.keys())
	
	counter = gen.OperationGenerator(
		a = map(zero_data, region_stream),
		b = count_stream,
		match = lambda a,b: b.start.shifted_forward(position) in a,
		a_is_passed = lambda a,b: b - a > abs(position), # this is a cautious maximum distance
		b_is_passed = lambda a,b: a - b > abs(position), # this is a cautious maximum distance
		operate = increment_count
	)
	
	if not quiet: bar = gen.ProgressBar(reference_lengths)
	for region in counter:
		if zeroes or region.data['count'] > 0:
			try:
				region_name = region.data['name']
			except KeyError:
				region_name = '%s:%i-%i%s' % (reference_names[region.reference_id], region.left_pos, region.right_pos, ('-' if region.is_reverse else '+'))
			out_tsv.write('%s\t%i\n' % (region_name, region.data['count']))
		if not quiet: bar.update(region)
	if not quiet: bar.finish()

if __name__ == '__main__':
	region_counts()

