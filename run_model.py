from tdc_models import TumorDeconvolutionModels
import argparse

parser = argparse.ArgumentParser(description='DREAM Tumor Deconvolution Challenge, team mitten_TDC19')

parser.add_argument('--input', type=str, default=None, help='input spreadsheet')
parser.add_argument('--model-level', type=str, default=None, help='model level, coarse or fine')
#parser.add_argument('--number-hvg', type=int, default=None, help='keep top marker genes with a high coefficient of variation, e.g. 15')
parser.add_argument('--signature-matrix', type=str, default=None, help='name of the signature matrix spreadsheet file')
#parser.add_argument('--model', type=str, default='median_random_samples', help='model name: arithmetic_mean, median_random_samples, linear_regression, elastic_net, arithmetic_mean_minus, geometric_mean, or log_2_geometric_mean')

args = parser.parse_args()

TumorDeconvolutionModels(model_level=args.model_level,
                         input_name=args.input,
                         model_name='median_random_samples',
                         use_gold_standard=False,
                         signature_matrix=args.signature_matrix).run()
