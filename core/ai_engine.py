# core/ai_engine.py
import streamlit as st
import torch
import re
import math
from transformers import AutoTokenizer, AutoModelForMaskedLM
from config import MODEL_NAME, KNOWN_MUTATIONS

class MutationPredictor:
    def __init__(self):
        self.tokenizer, self.model = self._load_resources()

    @staticmethod
    @st.cache_resource
    def _load_resources():
        """Cached loader for heavy model weights."""
        tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
        model = AutoModelForMaskedLM.from_pretrained(MODEL_NAME)
        return tokenizer, model

    def compute_score(self, seq: str, mutation: str) -> float:
        """Computes the log-likelihood ratio (Delta score)."""
        # Regex parsing
        m = re.match(r"([A-Z])(\d+)([A-Z])", mutation.upper())
        if not m:
            raise ValueError("Invalid mutation format. Use format like A23T.")
        
        wt, pos_str, mut = m.groups()
        pos = int(pos_str) - 1

        if pos < 0 or pos >= len(seq):
            raise ValueError(f"Position {pos+1} is out of bounds for sequence length {len(seq)}.")

        # Masking and Inference
        masked_seq = seq[:pos] + "<mask>" + seq[pos+1:]
        inputs = self.tokenizer(masked_seq, return_tensors="pt")
        
        with torch.no_grad():
            logits = self.model(**inputs).logits

        # Calculate probabilities
        mask_token_index = (inputs.input_ids == self.tokenizer.mask_token_id)[0].nonzero(as_tuple=True)[0]
        probs = torch.nn.functional.softmax(logits[0, mask_token_index, :], dim=-1)
        
        wt_id = self.tokenizer.convert_tokens_to_ids(wt)
        mut_id = self.tokenizer.convert_tokens_to_ids(mut)
        
        wt_prob = probs[0, wt_id].item()
        mut_prob = probs[0, mut_id].item()

        if wt_prob == 0:
            return 0.0
            
        return float(math.log(mut_prob / wt_prob))

    def analyze_impact(self, delta: float, mutation: str, gene: str) -> dict:
        """Combines AI score with Knowledge Base."""
        key = f"{gene}_{mutation.upper()}"
        
        # 1. Check Hardcoded Knowledge Base
        if key in KNOWN_MUTATIONS:
            return {
                "score": delta,
                "severity": KNOWN_MUTATIONS[key]['severity'],
                "insight": KNOWN_MUTATIONS[key]['summary'],
                "source": "Knowledge Base"
            }

        # 2. AI Prediction Heuristics
        if delta < -0.5:
            severity = "🔴 Likely Deleterious"
        elif delta < 0:
            severity = "🟠 Possibly Harmful"
        else:
            severity = "🟢 Likely Benign"

        return {
            "score": delta,
            "severity": severity,
            "insight": f"The AI model predicts a score of {delta:.3f}. Negative scores suggest the mutation is less probable than the wild type.",
            "source": "AI Model"
        }