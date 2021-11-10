#!/usr/bin/python
"""General purpose module ta handle blast file objects."""
import split_file as sp


class BlastFile:
    """Something"""
    def __init__(self, blast_name):
        self.name = blast_name
        self.blfi = open(blast_name)
        self.blast_ver = ''
        self.db = ''
        self.nr_queries = 0
        self.hit_percent = 0

        def read_pars(blast_file):
            no_hit_cnt = 0
            for lines in blast_file:
                if self.blast_ver == '':
                    self.blast_ver = lines.strip()
                if self.db == '' and lines.split(':')[0].strip() == 'Database':
                    self.db = lines.split(':')[1].strip()[:-1]
                if lines.split(' ')[0] == 'Query=':
                    self.nr_queries += 1
                if lines.strip() == '***** No hits found ******':
                    no_hit_cnt += 1
            self.hit_percent = round(100*(self.nr_queries - no_hit_cnt)
                                     / self.nr_queries, 1)
        read_pars(self.blfi)

    def print_par(self):
        print(self.nr_queries)
        print('Nr of queries: {}\nHit percentage: {} %'.format(
            self.nr_queries, self.hit_percent))
        print('Database: {}\nBlast version: {}'.format(self.db, self.blast_ver))

    def split(self, nr_splits):
        self.blfi.seek(0)
        sep = 'BLAST'
        return sp.split(self.blfi, sep, nr_splits)









