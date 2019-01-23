#!/usr/bin/python
"""General purpose module ta handle blast file objects."""


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
            self.hit_percent = round(100*(self.nr_queries - no_hit_cnt)/self.nr_queries, 1)
        read_pars(self.blfi)

    def print_par(self):
        print(self.nr_queries)
        print('Nr of queries: {}\nHit percentage: {} %'.format(self.nr_queries, self.hit_percent))
        print('Database: {}\nBlast version: {}'.format(self.db, self.blast_ver))

    def split(self, n):
        self.blfi.seek(0)
        bl_files = []
        if '.' in self.name:
            basename = self.name.split('.')[0]
        else:
            basename = self.name
        for count in range(n):
            bl_files.append(open(basename + '_' + str(count) + '.blast', 'w'))
        file_nr = 0
        query_counter = 1
        res_counter = self.nr_queries%n
        for line in self.blfi:
            try:
                if query_counter <= (self.nr_queries//n + (res_counter > 0)):
                    bl_files[file_nr].write(line)
                    if 'Query=' in line:
                        query_counter += 1
                elif self.blast_ver not in line:
                    bl_files[file_nr].write(line)
                else:
                    file_nr += 1
                    query_counter = 1
                    res_counter -= 1
                    bl_files[file_nr].write(line)
            except EOFError:
                for count in range(n):
                    bl_files[n].close()
        for count in range(n):
            bl_files.append(open(basename + '_' + str(count) + '.blast'))
        return bl_files








