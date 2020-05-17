import numpy as np

from GOTool import GeneOntology
from Utils import ColourClass, FancyApp


class BuildMatrices(FancyApp.FancyApp):

    def __init__(
            self,
            interpro,
            interpro_diff,
            interpro_cm,
            hmmer,
            hmmer_diff,
            hmmer_cm,
            s2f,
            goa,
            fasta,
            obo,
            domains=['cellular_component', 'biological_process',
                     'molecular_function'],
            verbose=True):

        super(BuildMatrices, self).__init__()
        self.colour = ColourClass.bcolors.BOLD_CYAN

        self.interpro = open(interpro, 'r')
        self.interpro_values = None

        self.interpro_diff = open(interpro_diff, 'r')
        self.interpro_diff_values = None

        if interpro_cm is not None:
            self.interpro_cm = open(interpro_cm, 'r')
            self.interpro_cm_values = None
        else:
            self.interpro_cm = None

        self.hmmer = open(hmmer, 'r')
        self.hmmer_values = None
        self.hmmer_diff = open(hmmer_diff, 'r')
        self.hmmer_diff_values = None

        if hmmer_cm is not None:
            self.hmmer_cm = open(hmmer_cm, 'r')
            self.hmmer_cm_values = None
        else:
            self.hmmer_cm = None

        self.s2f = open(s2f, 'r')
        self.s2f_values = None

        self.goa = open(goa, 'r')
        self.goa_values = None

        self.fasta = open(fasta, 'r')

        self.obo = open(obo, 'r')
        self.domains = domains

        self.go_terms = set()
        self.genes = set()

        self.__verbose__ = verbose

        evidence_codes = ['EXP', 'IDA', 'IPI', 'IMP',
                          'IGI', 'IEP', 'TAS', 'IC']
        self.gn = GeneOntology.GeneOntology(self.obo.name, verbose=True)
        self.gn.build_structure()
        self.gn.load_annotation_file(self.goa.name, 's2f', evidence_codes)
        self.gn.up_propagate_annotations('s2f')
        self.information_content = None

        self.infer_headers()
        self.infer_values()
        self.compute_information_content()

    def compute_information_content(self):
        for term, idx in self.go_terms_idx.items():
            self.information_content[idx] = self.gn.find_term(term)\
                                                .information_content('s2f')

    def inner_infer_headers(self, fp):
        for line in fp:
            parts = line.split()
            self.genes.add(parts[0])
            self.go_terms.add(parts[1])

    def infer_headers(self):
        self.tell('Inferring genes and terms from InterPro')
        self.inner_infer_headers(self.interpro)

        self.tell('Inferring genes and terms from HMMer')
        self.inner_infer_headers(self.hmmer)

        self.tell('Inferring genes and terms from InterPro + Diffusion')
        self.inner_infer_headers(self.interpro_diff)

        if self.interpro_cm is not None:
            self.tell('Inferring genes and terms from InterPro + CM Diffusion')
            self.inner_infer_headers(self.interpro_cm)

        self.tell('Inferring genes and terms from HMMer + Diffusion')
        self.inner_infer_headers(self.hmmer_diff)

        if self.hmmer_cm is not None:
            self.tell('Inferring genes and terms from HMMer + CM Diffusion')
            self.inner_infer_headers(self.hmmer_cm)

        self.tell('Inferring genes and terms from S2F')
        self.inner_infer_headers(self.s2f)

        self.tell('Inferring genes and terms from GOA')
        annotations = self.gn.get_annotations('s2f')
        self.genes |= set(annotations['Protein'].unique())
        self.go_terms |= set(annotations['GO ID'].unique())

        # transform the sets into ordered lists
        self.tell('Sorting terms...')
        self.go_terms = sorted(list(self.go_terms))

        self.tell('Indexing terms...')
        self.go_terms_idx = {}
        idx = 0
        for term in self.go_terms:
            self.go_terms_idx[term] = idx
            idx += 1

        self.tell('Sorting genes...')
        self.genes = sorted(list(self.genes))

        self.tell('Indexing genes...')
        self.genes_idx = {}
        idx = 0
        for gene in self.genes:
            self.genes_idx[gene] = idx
            idx += 1

        self.interpro_values = np.zeros((len(self.genes), len(self.go_terms)))
        self.hmmer_values = np.zeros((len(self.genes), len(self.go_terms)))
        self.interpro_diff_values = np.zeros((len(self.genes),
                                              len(self.go_terms)))
        self.hmmer_diff_values = np.zeros((len(self.genes),
                                           len(self.go_terms)))

        if self.interpro_cm is not None:
            self.interpro_cm_values = np.zeros((len(self.genes),
                                                len(self.go_terms)))
        if self.hmmer_cm is not None:
            self.hmmer_cm_values = np.zeros((len(self.genes),
                                             len(self.go_terms)))
        self.s2f_values = np.zeros((len(self.genes), len(self.go_terms)))
        self.goa_values = np.zeros((len(self.genes), len(self.go_terms)))
        self.information_content = np.zeros(len(self.go_terms))

    def inner_infer_values_goa(self, value_array):
        annotations = self.gn.get_annotations('s2f')
        for i, row in annotations.iterrows():
            value_array[self.genes_idx[row['Protein']],
                        self.go_terms_idx[row['GO ID']]] = 1.0

    def inner_infer_values(self, fp, value_array):
        fp.seek(0)
        for line in fp:
            parts = line.split()
            if float(parts[2]) > 0:
                value_array[self.genes_idx[parts[0]],
                            self.go_terms_idx[parts[1]]] = float(parts[2])

    def infer_values(self):
        self.tell('Inferring values from InterPro')
        self.inner_infer_values(self.interpro, self.interpro_values)

        self.tell('Inferring values from HMMer')
        self.inner_infer_values(self.hmmer, self.hmmer_values)

        self.tell('Inferring values from InterPro + Diffusion')
        self.inner_infer_values(self.interpro_diff, self.interpro_diff_values)

        if self.interpro_cm is not None:
            self.tell('Inferring values from InterPro + CM Diffusion')
            self.inner_infer_values(self.interpro_cm, self.interpro_cm_values)

        self.tell('Inferring values from HMMer + Diffusion')
        self.inner_infer_values(self.hmmer_diff, self.hmmer_diff_values)

        if self.hmmer_cm is not None:
            self.tell('Inferring values from HMMer + CM Diffusion')
            self.inner_infer_values(self.hmmer_cm, self.hmmer_cm_values)

        self.tell('Inferring values from S2F')
        self.inner_infer_values(self.s2f, self.s2f_values)

        self.tell('Inferring values from GOA')
        self.inner_infer_values_goa(self.goa_values)
