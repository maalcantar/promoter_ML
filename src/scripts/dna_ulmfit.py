import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
# this to set the plot styles 
sns.set_style("darkgrid", {"axes.facecolor": ".9"})
sns.set_context('paper', font_scale = 1.5)
from matplotlib.colors import LogNorm

from Bio import Seq, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation
import networkx as nx
from itertools import islice
import re

import torch
from torch.utils.data.sampler import WeightedRandomSampler
from fastai import *
from fastai.text import *
from fastai.text.interpret import *
from fastai.callbacks import *
__all__ = ['OverSamplingCallback']

from ranger import Ranger
optar = partial(Ranger)

# Ignore  the warnings
import warnings
warnings.filterwarnings('always')
warnings.filterwarnings('ignore')

# import sklearn modules 
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_curve, auc, confusion_matrix
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.feature_extraction.text import CountVectorizer, TfidfTransformer, TfidfVectorizer
from sklearn.utils.multiclass import unique_labels

label_encoder = LabelEncoder()
kfold = StratifiedKFold(n_splits = 5, shuffle = True)
# implementing a bag of words model where we weight each word by its document frequency 
tfidf_vec = TfidfVectorizer(ngram_range = (1,3), lowercase = True) # these are the ngram parameters to sweep over 

import tensorflow as tf
tf.config.optimizer.set_jit(True)

from gensim.utils import *
from multiprocessing import cpu_count
cores = multiprocessing.cpu_count() # Count the number of cores in a computer

# define functions to k-merize DNA and produce a 'sentence' of k-mer tokens 
def dna_kmerize(seq, ngram_size, ngram_stride):
    if ngram_size == 1:
        toks = list(seq) # for character level encoding of a dna sequence
    else:
        toks = [seq[i:i+ngram_size] for i in range(0, len(seq), ngram_stride) if len(seq[i:i+ngram_size]) == ngram_size]
    if len(toks[-1]) < ngram_size: # throw out tokens that don't match the length of the kmer (as opposed to a pad token, which serves the same purpose)
        toks = toks[:-1]
    return toks

# stitch the tokens together into a sentence 
def seq_sentence(seq, ngram_size, ngram_stride):
    
    toks = dna_kmerize(seq, ngram_size, ngram_stride)
    seq_sentence = ' '.join(toks)
    return seq_sentence


# this function takes a list of strings and prepares them into a databunch with identifying tokens 
def _join_dna_texts(texts:Collection[str],
    mark_fields:bool = False,
    include_bos:bool = False,
    include_eos:bool = False):
  
    if not isinstance(texts, np.ndarray): 
      texts = np.array(texts)
    if is1d(texts): 
      texts = texts[:,None]
    df = pd.DataFrame({i:texts[:,i] for i in range(texts.shape[1])})

    bos_tok = f'{BOS} ' if include_bos else ''
    text_col = f'{bos_tok}{FLD} {1} ' + df[0].astype(str) if mark_fields else f'{bos_tok}' + df[0].astype(str)

    for i in range(1,len(df.columns)):
        text_col += (f' {FLD} {i+1} ' if mark_fields else ' ') + df[i].astype(str)

    if include_eos: 
      text_col = text_col + f' {EOS}'

    return text_col.values


# we call these GenomicTokenizeProcessor and GenomicNumericalizeProcessor respectively 
# this creates a vocabulary with associated dictionary building off the fastai Vocab class
# so we don't need to create the other associated functions 
class GenomicVocab(Vocab):
    "Contain the correspondence between numbers and tokens and numericalizer."
    def __init__(self, 
        itos:Collection[str]):

        self.itos = itos
        self.stoi = collections.defaultdict(int,{v:k for k,v in enumerate(self.itos)})

    @classmethod
    def create(cls, 
        tokens:Tokens, 
        max_vocab:int, 
        min_freq:int) -> 'Vocab':
        "Create a vocabulary from a set of `tokens`. itos is a list that contains the id to kmer correspondence."

        freq = Counter(p for o in tokens for p in o)
        itos = [o for o,c in freq.most_common(max_vocab) if c >= min_freq]
        itos.insert(0, 'pad') # all pad tokens have the same integer id 0
        return cls(itos)

class GenomicTokenizeProcessor(PreProcessor):
    "PreProcessor that tokenizes the texts in ds."

    def __init__(self, 
        ds:ItemList = None, 
        tokenizer:Tokenizer = None, 
        chunksize:int = 10000, 
        mark_fields:bool = False, 
        include_bos:bool = False, 
        include_eos:bool = False):
      
        self.tokenizer, self.chunksize, self.mark_fields = ifnone(tokenizer, Tokenizer()), chunksize, mark_fields
        self.include_bos, self.include_eos = include_bos, include_eos

    def process_one(self, item):
        return self.tokenizer._process_all_1(_join_dna_texts([item], self.mark_fields, self.include_bos, self.include_eos))[0]

    def process(self, ds):
        ds.items = _join_dna_texts(ds.items, self.mark_fields, self.include_bos, self.include_eos)

        tokens = []
        for i in progress_bar(range(0,len(ds),self.chunksize), leave = False):
            tokens += self.tokenizer.process_all(ds.items[i:i+self.chunksize])
        ds.items = tokens

class GenomicNumericalizeProcessor(PreProcessor):
    "PreProcessor that numericalizes the tokens in ds."

    def __init__(self, 
        ds:ItemList = None, 
        vocab:Vocab = None, 
        max_vocab:int = 60000, 
        min_freq:int = 2):

        vocab = ifnone(vocab, ds.vocab if ds is not None else None)
        self.vocab, self.max_vocab, self.min_freq = vocab, max_vocab, min_freq

    def process_one(self, item): 
      return np.array(self.vocab.numericalize(item), dtype = np.int64)

    def process(self, ds):
        if self.vocab is None: 
            self.vocab = GenomicVocab.create(ds.items, self.max_vocab, self.min_freq)
        ds.vocab = self.vocab
        super().process(ds)

# this function is essentially a pipeline that executes the tokenization and numericalization for the databunch routine 
# we need a custom tokenizer and enumerator for tokenizing and integer indexing genomes.
def _get_genomic_processor(tokenizer:Tokenizer = None,
    vocab:Vocab = None, 
    chunksize:int = 10000, 
    max_vocab:int = 60000,
    min_freq:int = 2, 
    mark_fields:bool = False, 
    include_bos:bool = False, 
    include_eos:bool = False):
  
    return [
        GenomicTokenizeProcessor(tokenizer = tokenizer, chunksize = chunksize, mark_fields = mark_fields, include_bos = include_bos, include_eos = include_eos),
        GenomicNumericalizeProcessor(vocab = vocab, max_vocab = max_vocab, min_freq = min_freq)
            ]

# now we can put this all together for the textclassification and textLM databunch routines from df
class GenomicTextClasDataBunch(TextClasDataBunch):
    "Create a TextDataBunch suitable for training an RNN classifier."
    @classmethod
    def from_df(cls, 
        path:PathOrStr, 
        train_df:DataFrame, 
        valid_df:DataFrame, 
        test_df:Optional[DataFrame] = None,
        tokenizer:Tokenizer = None, 
        vocab:Vocab = None, 
        classes:Collection[str] = None, 
        text_cols:IntsOrStrs = 1,
        label_cols:IntsOrStrs = 0, 
        label_delim:str = None, 
        chunksize:int = 10000, 
        max_vocab:int = 60000,
        min_freq:int = 2, # take all terms 
        mark_fields:bool = False, 
        include_bos:bool = False,
        include_eos:bool = False, 
        bs:int = 32, # batch_size 
        val_bs:int = None, 
        pad_idx = 0, # keep this the same as this is the value assigned to the pad token 
        pad_first = True, 
        device:torch.device = None, # use a CPU for preparing databunch routine
        no_check:bool = False, 
        backwards:bool = False, **kwargs) -> DataBunch: # by default backwards is false
        
        "Create a TextDataBunch from DataFrames. `kwargs` are passed to the dataloader creation. The first part is copied from TextDataBunch from_df method"
        processor = _get_genomic_processor(
            tokenizer = tokenizer, 
            vocab = vocab, 
            chunksize = chunksize, 
            max_vocab = max_vocab, 
            min_freq = min_freq, 
            mark_fields = mark_fields, 
            include_bos = include_bos, 
            include_eos = include_eos
            )

        if classes is None and is_listy(label_cols) and len(label_cols) > 1: 
            classes = label_cols

        src = ItemLists(path, 
            TextList.from_df(train_df, path, cols = text_cols, processor = processor),
            TextList.from_df(valid_df, path, cols = text_cols, processor = processor)
            )
        
        if label_delim is not None: 
          src = src.label_from_df(cols = label_cols, classes = classes, label_delim = label_delim)
        else: 
          src = src.label_from_df(cols = label_cols, classes = classes)

        if test_df is not None: 
          src.add_test(TextList.from_df(test_df, path, cols = text_cols))

        dl = src.databunch(**kwargs) # this is now the dataloader which will be fed into the databunch routine 

        # "Function that transform the `datasets` in a `DataBunch` for classification. Passes `**dl_kwargs` on to `DataLoader()`"
        datasets = cls._init_ds(dl.train_ds, dl.valid_ds, dl.test_ds)
        val_bs = ifnone(val_bs, bs)
        collate_fn = partial(pad_collate, pad_idx = pad_idx, pad_first = pad_first, backwards = backwards)

        train_sampler = SortishSampler(datasets[0].x, key=lambda t: len(datasets[0][t][0].data), bs = bs) # mixes up the text-order 
        train_dl = DataLoader(datasets[0], batch_size = bs, sampler = train_sampler, drop_last = True, **kwargs)
        dataloaders = [train_dl]

        for ds in datasets[1:]:
            lengths = [len(t) for t in ds.x.items]
            sampler = SortSampler(ds.x, key = lengths.__getitem__)
            dataloaders.append(DataLoader(ds, batch_size = val_bs, sampler = sampler, **kwargs))

        return cls(*dataloaders, path = path, collate_fn = collate_fn, no_check = no_check)

class GenomicTextLMDataBunch(TextLMDataBunch):
    "Create a `TextDataBunch` suitable for training a language model. Combine TextDataBunch and TextLMDatabunch" 
    @classmethod
    def from_df(cls, 
        path:PathOrStr, 
        train_df:DataFrame, 
        valid_df:DataFrame, 
        test_df:Optional[DataFrame] = None,
        tokenizer:Tokenizer = None, 
        vocab:Vocab = None, 
        classes:Collection[str] = None, 
        text_cols:IntsOrStrs = 1,
        label_cols:IntsOrStrs = 0, 
        label_delim:str = None, 
        chunksize:int = 10000, 
        max_vocab:int = 60000,
        min_freq:int = 2, 
        mark_fields:bool = False, 
        include_bos:bool = False, 
        include_eos:bool = False, 
        no_check:bool = False, 
        bs = 64, 
        val_bs:int = None, 
        num_workers:int = 0, 
        device:torch.device = None, 
        collate_fn:Callable = data_collate, 
        bptt:int = 30, # back-prop through time 
        backwards:bool = False, **kwargs) -> DataBunch:
        
        processor = _get_genomic_processor(
            tokenizer = tokenizer, 
            vocab = vocab, 
            chunksize = chunksize, 
            max_vocab = max_vocab,
            min_freq = min_freq, 
            mark_fields = mark_fields, 
            include_bos = include_bos, 
            include_eos = include_eos
            )
        
        if classes is None and is_listy(label_cols) and len(label_cols) > 1: 
            classes = label_cols

        src = ItemLists(path, 
            TextList.from_df(train_df, path, cols = text_cols, processor = processor),
            TextList.from_df(valid_df, path, cols = text_cols, processor = processor)
            )

        src = src.label_for_lm()

        if test_df is not None: 
            src.add_test(TextList.from_df(test_df, path, cols=text_cols))
        
        dl = src.databunch(**kwargs)

        datasets = cls._init_ds(dl.train_ds, dl.valid_ds, dl.test_ds)
        val_bs = ifnone(val_bs, bs)
        datasets = [LanguageModelPreLoader(ds, shuffle=(i==0), bs=(bs if i==0 else val_bs), bptt = bptt, backwards = backwards)
                    for i,ds in enumerate(datasets)]
        val_bs = bs
        dls = [DataLoader(d, b, shuffle = False, **kwargs) for d,b in zip(datasets, (bs, val_bs, val_bs, val_bs)) if d is not None]

        return cls(*dls, path = path, collate_fn = collate_fn, no_check = no_check)

class OverSamplingCallback(LearnerCallback):
    def __init__(self,
        learn:Learner,
        weights:torch.Tensor=None):

        super().__init__(learn)
        self.weights = weights

    def on_train_begin(self, **kwargs):
        ds, dl = self.data.train_ds, self.data.train_dl
        self.labels = ds.y.items

        assert np.issubdtype(self.labels.dtype, np.integer), "Can only oversample integer values"

        _, self.label_counts = np.unique(self.labels, return_counts = True)

        if self.weights is None: 
            self.weights = torch.DoubleTensor((1/self.label_counts)[self.labels])

        self.total_len_oversample = int(self.data.c*np.max(self.label_counts))
        sampler = WeightedRandomSampler(self.weights, self.total_len_oversample)
        self.data.train_dl = dl.new(shuffle = False, sampler = sampler)

# lets define the functions that make the language models 
_model_meta = {AWD_LSTM: {'hid_name':'emb_sz', 'url':URLs.WT103_FWD, 'url_bwd':URLs.WT103_BWD,
                          'config_lm':awd_lstm_lm_config, 'split_lm': awd_lstm_lm_split,
                          'config_clas':awd_lstm_clas_config, 'split_clas': awd_lstm_clas_split},
               Transformer: {'hid_name':'d_model', 'url':URLs.OPENAI_TRANSFORMER,
                             'config_lm':tfmer_lm_config, 'split_lm': tfmer_lm_split,
                             'config_clas':tfmer_clas_config, 'split_clas': tfmer_clas_split},
               TransformerXL: {'hid_name':'d_model',
                              'config_lm':tfmerXL_lm_config, 'split_lm': tfmerXL_lm_split,
                              'config_clas':tfmerXL_clas_config, 'split_clas': tfmerXL_clas_split}}

def get_dnalanguage_model(
    arch:Callable, 
    vocab_sz:int, 
    config:dict = None, 
    drop_mult:float = 1.):
    "Create a language model from `arch` and its `config`, maybe `pretrained`."

    _model_meta = {AWD_LSTM: {'hid_name':'emb_sz', 'url':URLs.WT103_FWD, 'url_bwd':URLs.WT103_BWD,
                          'config_lm':awd_lstm_lm_config, 'split_lm': awd_lstm_lm_split,
                          'config_clas':awd_lstm_clas_config, 'split_clas': awd_lstm_clas_split},
               Transformer: {'hid_name':'d_model', 'url':URLs.OPENAI_TRANSFORMER,
                             'config_lm':tfmer_lm_config, 'split_lm': tfmer_lm_split,
                             'config_clas':tfmer_clas_config, 'split_clas': tfmer_clas_split},
               TransformerXL: {'hid_name':'d_model',
                              'config_lm':tfmerXL_lm_config, 'split_lm': tfmerXL_lm_split,
                              'config_clas':tfmerXL_clas_config, 'split_clas': tfmerXL_clas_split}}

    meta = _model_meta[arch]
    config = ifnone(config, meta['config_lm']).copy()

    for k in config.keys():
        if k.endswith('_p'): config[k] *= drop_mult # multiply all dropout with multiplier 

    tie_weights, output_p, out_bias = map(config.pop, ['tie_weights', 'output_p', 'out_bias'])
    init = config.pop('init') if 'init' in config else None
    encoder = arch(vocab_sz, **config)
    enc = encoder.encoder if tie_weights else None
    decoder = LinearDecoder(vocab_sz, config[meta['hid_name']], output_p, tie_encoder = enc, bias = out_bias)
    model = SequentialRNN(encoder, decoder)
    return model if init is None else model.apply(init)

def dnalanguage_model_learner(
    data:DataBunch, 
    arch, 
    config:dict = None, 
    drop_mult:float = 1., 
    pretrained:bool = False,
    pretrained_fnames:OptStrTuple = None, **learn_kwargs) -> 'LanguageLearner':
    "Create a `Learner` with a language model from `data` and `arch`."
    _model_meta = {AWD_LSTM: {'hid_name':'emb_sz', 'url':URLs.WT103_FWD, 'url_bwd':URLs.WT103_BWD,
                          'config_lm':awd_lstm_lm_config, 'split_lm': awd_lstm_lm_split,
                          'config_clas':awd_lstm_clas_config, 'split_clas': awd_lstm_clas_split},
               Transformer: {'hid_name':'d_model', 'url':URLs.OPENAI_TRANSFORMER,
                             'config_lm':tfmer_lm_config, 'split_lm': tfmer_lm_split,
                             'config_clas':tfmer_clas_config, 'split_clas': tfmer_clas_split},
               TransformerXL: {'hid_name':'d_model',
                              'config_lm':tfmerXL_lm_config, 'split_lm': tfmerXL_lm_split,
                              'config_clas':tfmerXL_clas_config, 'split_clas': tfmerXL_clas_split}}

    model = get_dnalanguage_model(
        arch, 
        len(data.vocab.itos), 
        config = config, 
        drop_mult = drop_mult
        )
    meta = _model_meta[arch]
    learn = LanguageLearner(data, model, split_func=meta['split_lm'], 
                            callback_fns = [partial(EarlyStoppingCallback, monitor = 'accuracy', min_delta = 0.001, patience = 3)], **learn_kwargs)

    url = 'url_bwd' if data.backwards else 'url'
    if pretrained or pretrained_fnames:
        if pretrained_fnames is not None:
            fnames = [learn.path/learn.model_dir/f'{fn}.{ext}' for fn,ext in zip(pretrained_fnames, ['pth', 'pkl'])]
        else:
            if url not in meta:
                warn("There are no pretrained weights for that architecture yet!")
                return learn
            model_path = untar_data(meta[url] , data=False)
            fnames = [list(model_path.glob(f'*.{ext}'))[0] for ext in ['pth', 'pkl']]
        learn = learn.load_pretrained(*fnames)
        learn.freeze()

    return learn

def get_dna_classifier(
    arch:Callable, 
    vocab_sz:int, 
    n_class:int, 
    bptt:int = 30, 
    max_len:int = 20*30, 
    config:dict = None,
    drop_mult:float = 1., 
    lin_ftrs:Collection[int] = None, 
    ps:Collection[float] = None,
    pad_idx:int=1) -> nn.Module:
    "Create a text classifier from `arch` and its `config`, maybe `pretrained`."
    _model_meta = {AWD_LSTM: {'hid_name':'emb_sz', 'url':URLs.WT103_FWD, 'url_bwd':URLs.WT103_BWD,
                          'config_lm':awd_lstm_lm_config, 'split_lm': awd_lstm_lm_split,
                          'config_clas':awd_lstm_clas_config, 'split_clas': awd_lstm_clas_split},
               Transformer: {'hid_name':'d_model', 'url':URLs.OPENAI_TRANSFORMER,
                             'config_lm':tfmer_lm_config, 'split_lm': tfmer_lm_split,
                             'config_clas':tfmer_clas_config, 'split_clas': tfmer_clas_split},
               TransformerXL: {'hid_name':'d_model',
                              'config_lm':tfmerXL_lm_config, 'split_lm': tfmerXL_lm_split,
                              'config_clas':tfmerXL_clas_config, 'split_clas': tfmerXL_clas_split}}

    meta = _model_meta[arch]
    config = ifnone(config, meta['config_clas']).copy()

    for k in config.keys():
        if k.endswith('_p'): config[k] *= drop_mult

    if lin_ftrs is None: 
        lin_ftrs = [50]
    if ps is None:  
        ps = [0.1]*len(lin_ftrs)

    layers = [config[meta['hid_name']] * 3] + lin_ftrs + [n_class]
    ps = [config.pop('output_p')] + ps
    init = config.pop('init') if 'init' in config else None

    encoder = MultiBatchEncoder(bptt, max_len, arch(vocab_sz, **config), pad_idx = pad_idx)
    model = SequentialRNN(encoder, PoolingLinearClassifier(layers, ps))
    return model if init is None else model.apply(init)

def dna_classifier_learner(
    data:DataBunch, 
    arch:Callable, 
    bptt:int = 30, 
    max_len:int = 20*30, 
    config:dict = None,
    pretrained:bool = False, 
    drop_mult:float = 1., 
    lin_ftrs:Collection[int] = None,
    ps:Collection[float] = None, **learn_kwargs) -> 'TextClassifierLearner':
    "Create a `Learner` with a text classifier from `data` and `arch`."

    _model_meta = {AWD_LSTM: {'hid_name':'emb_sz', 'url':URLs.WT103_FWD, 'url_bwd':URLs.WT103_BWD,
                          'config_lm':awd_lstm_lm_config, 'split_lm': awd_lstm_lm_split,
                          'config_clas':awd_lstm_clas_config, 'split_clas': awd_lstm_clas_split},
               Transformer: {'hid_name':'d_model', 'url':URLs.OPENAI_TRANSFORMER,
                             'config_lm':tfmer_lm_config, 'split_lm': tfmer_lm_split,
                             'config_clas':tfmer_clas_config, 'split_clas': tfmer_clas_split},
               TransformerXL: {'hid_name':'d_model',
                              'config_lm':tfmerXL_lm_config, 'split_lm': tfmerXL_lm_split,
                              'config_clas':tfmerXL_clas_config, 'split_clas': tfmerXL_clas_split}}

    model = get_dna_classifier(
                arch, 
                len(data.vocab.itos), 
                data.c, 
                bptt = bptt, 
                max_len = max_len,
                config = config, 
                drop_mult = drop_mult, 
                lin_ftrs = lin_ftrs, 
                ps = ps)

    meta = _model_meta[arch]
    learn = RNNLearner(
                data, 
                model, 
                split_func = meta['split_clas'],
                callback_fns = [partial(EarlyStoppingCallback, monitor = 'accuracy', min_delta = 0.01, patience = 3)],
                opt_func = optar, 
                bn_wd = False, 
                true_wd = True,
                **learn_kwargs)
    #callback_fns=[partial(OverSamplingCallback)]
    if pretrained:
        if 'url' not in meta:
            warn("There are no pretrained weights for that architecture yet!")
            return learn
        model_path = untar_data(meta['url'], data=False)
        fnames = [list(model_path.glob(f'*.{ext}'))[0] for ext in ['pth', 'pkl']]
        learn = learn.load_pretrained(*fnames, strict=False)
        learn.freeze()

    return learn

