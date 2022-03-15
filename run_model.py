model_level = 'fine' # 'coarse' or 'fine'
new_round = True

from tdc_models import TumorDeconvolutionModels

models = {1: 'arithmetic_mean',
          2: 'median_random_samples',
          3: 'linear_regression',
          4: 'elastic_net',
          5: 'arithmetic_mean_minus',
          6: 'geometric_mean',
          7: 'log_2_geometric_mean'}


model = models[1]
sig = 'signature_matrix_combo_%s.xlsx' % (model_level)
print('Model name: %s, Sub-challenge: %s' % (model, model_level))
print('Using signature matrix: %s\n' % (sig))

if new_round:
    # Validation round
    num2keep = 15
    TumorDeconvolutionModels(model_level=model_level, input_name='input.csv', model_name=model, use_gold_standard=False, signature_matrix=sig).run(num2keep)
else:
    # keep_per_cell = [i for i in range(5, 30, 3)]
    keep_per_cell = [15]
    for num2keep in keep_per_cell:
        print("markers per celltype:", num2keep, "\n")
        for model_level in ["fine"]:
            sig = 'signature_matrix_combo_%s.xlsx' % (model_level)

            # Round 1
            TumorDeconvolutionModels(model_level=model_level, input_name='input_r1.csv', model_name=model, round=1, signature_matrix=sig).run(num2keep)

            # Round 2
            TumorDeconvolutionModels(model_level=model_level, input_name='input_r2.csv', model_name=model, round=2, signature_matrix=sig).run(num2keep)

            # Round 3
            TumorDeconvolutionModels(model_level=model_level, input_name='input_r3_%s.csv' % (model_level), model_name=model, round=3, signature_matrix=sig).run(num2keep)
