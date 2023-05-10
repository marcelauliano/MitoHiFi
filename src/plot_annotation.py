import sys
import numpy as np
import PIL
from dna_features_viewer import BiopythonTranslator
from BCBio import GFF

class MyCustomTranslator(BiopythonTranslator):


    def compute_filtered_features(self, features):
        return [feature for feature in features if feature.type != "gene"] 


class MyCustomTranslatorMitos(BiopythonTranslator):

    def compute_filtered_features(self, features):
        return [feature for feature in features if feature.type == "gene"] 

    def compute_feature_label(self, feature):
        if feature.type == "gene":
            return feature.id.replace("gene_", "")

def plot_annotation(in_gb, out_gb=None):
    
    if not out_gb:
        out_gb = f"{in_gb}".replace("gb", "png")

    graphic_record = MyCustomTranslator().translate_record(in_gb)
    ax, _ = graphic_record.plot(figure_width=10, strand_in_label_threshold=7)
    ax.set_title(out_gb)
    ax.figure.savefig(out_gb)


def plot_annotation_mitos(in_gff, out_gff=None):
    
    if not out_gff:
        out_gff = f"{in_gff}".replace("gff", "png")
    
    graphic_record = MyCustomTranslatorMitos().translate_record(in_gff)
    ax, _ = graphic_record.plot(figure_width=10, strand_in_label_threshold=7)
    ax.set_title(out_gff)
    ax.figure.savefig(out_gff)

def merge_images(img_list, out_file):

    imgs = [ PIL.Image.open(i) for i in img_list ]
    imgs_comb = np.vstack([np.asarray(i) for i in imgs])
    imgs_comb = PIL.Image.fromarray(imgs_comb)
    imgs_comb.save(out_file)




if __name__ == "__main__":
    
    out_gb = None

    num_args = len(sys. argv) - 1
    if num_args == 1:
        in_gb = sys.argv[1]
    elif num_args == 2:
        in_gb = sys.argv[1]
        out_gb = sys.argv[2]
    else:
        print("Error. Give at least one argument (input genbank filename) or at most two (input genbank filename and output plot name")

    plot_annotation(in_gb, out_gb)
