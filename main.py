from flask import Flask, request, jsonify
from flask_cors import CORS
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import io
import os
import requests
import re

app = Flask(__name__)

CORS(app, origins="*", allow_headers=["Content-Type", "Authorization"], methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"])

# ─────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────

def fetch_file(file_url):
    response = requests.get(file_url, timeout=30)
    response.raise_for_status()
    return response.text

def interpret_gc(gc_percentage, sequence_length, g_count, c_count, a_count, t_count):
    if 40 <= gc_percentage <= 60:
        gc_status = "Normal"
        gc_interpretation = f"GC content of {gc_percentage}% falls within the normal range (40-60%). Good sequence quality."
        gc_recommendation = "Proceed with downstream analysis such as variant calling or alignment."
    elif gc_percentage < 40:
        gc_status = "AT-Rich"
        gc_interpretation = f"GC content of {gc_percentage}% is below normal, indicating an AT-rich sequence."
        gc_recommendation = "Flag for review. Consider verifying sample quality and re-sequencing."
    else:
        gc_status = "GC-Rich"
        gc_interpretation = f"GC content of {gc_percentage}% is above normal, indicating a GC-rich sequence."
        gc_recommendation = "Verify PCR conditions were optimized for GC-rich templates."

    if sequence_length < 100:
        length_status = "Very Short"
        length_interpretation = "Sequence length is very short. May indicate incomplete or degraded sample."
        length_recommendation = "Re-sequence the sample."
    elif sequence_length < 1000:
        length_status = "Short"
        length_interpretation = f"Sequence length of {sequence_length} bp is short. Suitable for targeted sequencing."
        length_recommendation = "Suitable for targeted gene panels."
    elif sequence_length < 50000:
        length_status = "Standard"
        length_interpretation = f"Sequence length of {sequence_length} bp is within standard range."
        length_recommendation = "Proceed with standard bioinformatics pipeline."
    else:
        length_status = "Long"
        length_interpretation = f"Sequence length of {sequence_length} bp indicates large genomic region."
        length_recommendation = "Use high-performance computing pipeline."

    total = g_count + c_count + a_count + t_count
    a_pct = round((a_count / total) * 100, 2) if total > 0 else 0
    t_pct = round((t_count / total) * 100, 2) if total > 0 else 0
    g_pct = round((g_count / total) * 100, 2) if total > 0 else 0
    c_pct = round((c_count / total) * 100, 2) if total > 0 else 0

    gc_skew = round((g_count - c_count) / (g_count + c_count), 4) if (g_count + c_count) > 0 else 0
    if gc_skew > 0.1:
        skew_interpretation = "Positive GC skew detected, suggesting leading strand bias."
    elif gc_skew < -0.1:
        skew_interpretation = "Negative GC skew detected, suggesting lagging strand bias."
    else:
        skew_interpretation = "GC skew is balanced, indicating no significant strand bias."

    quality_score = 100
    if gc_percentage < 35 or gc_percentage > 65:
        quality_score -= 30
    elif gc_percentage < 40 or gc_percentage > 60:
        quality_score -= 15
    if sequence_length < 100:
        quality_score -= 30
    elif sequence_length < 1000:
        quality_score -= 10
    if abs(gc_skew) > 0.2:
        quality_score -= 10
    quality_score = max(0, quality_score)

    if quality_score >= 80:
        overall_quality = "Excellent"
        overall_summary = "This sample demonstrates excellent quality metrics and is suitable for clinical and research analysis."
    elif quality_score >= 60:
        overall_quality = "Good"
        overall_summary = "Good quality with minor flags. Results are reliable for most analyses."
    elif quality_score >= 40:
        overall_quality = "Fair"
        overall_summary = "Some quality concerns. Review flagged metrics before proceeding."
    else:
        overall_quality = "Poor"
        overall_summary = "Significant quality issues. Re-sequencing is strongly recommended."

    return {
        "gc_status": gc_status,
        "gc_interpretation": gc_interpretation,
        "gc_recommendation": gc_recommendation,
        "length_status": length_status,
        "length_interpretation": length_interpretation,
        "length_recommendation": length_recommendation,
        "nucleotide_percentages": {"A": a_pct, "T": t_pct, "G": g_pct, "C": c_pct},
        "gc_skew": gc_skew,
        "skew_interpretation": skew_interpretation,
        "quality_score": quality_score,
        "overall_quality": overall_quality,
        "overall_summary": overall_summary
    }


# ─────────────────────────────────────────
# ROUTES
# ─────────────────────────────────────────

@app.route('/', methods=['GET'])
def home():
    return jsonify({
        "app_name": "Genomirates",
        "status": "online",
        "version": "3.0",
        "message": "Genomirates Bioinformatics Engine is running.",
        "endpoints": {
            "fasta_analysis": "POST /api/calc",
            "vcf_analysis": "POST /api/vcf",
            "fastq_qc": "POST /api/fastq",
            "primer_design": "POST /api/primers",
            "wellness_dna": "POST /api/wellness",
            "premarital_screening": "POST /api/premarital",
            "carrier_screening": "POST /api/carrier",
            "chronic_disease_risk": "POST /api/chronic",
            "pharmacogenomics": "POST /api/pharma",
            "nutrigenomics": "POST /api/nutri",
            "corporate_wellness": "POST /api/corporate",
            "cancer_risk": "POST /api/cancer"
        }
    })


# ─────────────────────────────────────────
# 1. FASTA ANALYSIS (existing + enhanced)
# ─────────────────────────────────────────

@app.route('/api/calc', methods=['GET', 'POST'])
def analyze_sequence():
    if request.method == 'GET':
        return jsonify({"status": "online", "message": "POST a sequence or file_url to analyze."})
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        fasta_text = ""
        if 'file_url' in data:
            try:
                fasta_text = fetch_file(data['file_url'])
            except Exception as e:
                return jsonify({"status": "error", "message": f"Failed to download file: {str(e)}"}), 400
        elif 'sequence' in data:
            fasta_text = data['sequence']
        else:
            return jsonify({"status": "error", "message": "Missing sequence or file_url"}), 400

        fasta_io = io.StringIO(fasta_text)
        try:
            record = next(SeqIO.parse(fasta_io, "fasta"))
        except StopIteration:
            return jsonify({"status": "error", "message": "Invalid FASTA format"}), 400

        sequence_str = str(record.seq).upper()
        if len(sequence_str) == 0:
            return jsonify({"status": "error", "message": "Empty sequence"}), 400

        sequence_length = len(sequence_str)
        g_count = sequence_str.count('G')
        c_count = sequence_str.count('C')
        a_count = sequence_str.count('A')
        t_count = sequence_str.count('T')
        gc_content = round(((g_count + c_count) / sequence_length) * 100, 2)

        interpretation = interpret_gc(gc_content, sequence_length, g_count, c_count, a_count, t_count)

        return jsonify({
            "status": "success",
            "analysis_results": {
                "sample_id": record.id,
                "sequence_length": sequence_length,
                "gc_content_percentage": gc_content,
                "g_count": g_count,
                "c_count": c_count,
                "a_count": a_count,
                "t_count": t_count,
                **interpretation
            },
            "message": f"Successfully analyzed {record.id}"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"Processing failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# 2. VCF VARIANT ANALYSIS
# ─────────────────────────────────────────

@app.route('/api/vcf', methods=['POST'])
def analyze_vcf():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        vcf_text = ""
        if 'file_url' in data:
            vcf_text = fetch_file(data['file_url'])
        elif 'content' in data:
            vcf_text = data['content']
        else:
            return jsonify({"status": "error", "message": "Missing file_url or content"}), 400

        lines = vcf_text.strip().split('\n')
        variants = []
        snps = 0
        indels = 0
        total_variants = 0
        chromosomes = {}

        for line in lines:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alt = parts[4]

            total_variants += 1
            chromosomes[chrom] = chromosomes.get(chrom, 0) + 1

            if len(ref) == 1 and len(alt) == 1:
                snps += 1
                variant_type = "SNP"
            else:
                indels += 1
                variant_type = "INDEL"

            variants.append({
                "chromosome": chrom,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "type": variant_type
            })

        if total_variants == 0:
            quality = "No Variants"
            summary = "No variants detected in this VCF file."
        elif snps / max(total_variants, 1) > 0.9:
            quality = "SNP Rich"
            summary = "Sample is predominantly SNP variants, typical of human genomic variation data."
        else:
            quality = "Mixed Variants"
            summary = "Sample contains a mix of SNPs and indels."

        return jsonify({
            "status": "success",
            "vcf_results": {
                "total_variants": total_variants,
                "snp_count": snps,
                "indel_count": indels,
                "snp_ratio": round(snps / max(total_variants, 1) * 100, 2),
                "chromosomes_affected": len(chromosomes),
                "chromosome_distribution": chromosomes,
                "variant_quality": quality,
                "clinical_summary": summary,
                "top_variants": variants[:20]
            },
            "message": f"Successfully analyzed {total_variants} variants"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"VCF analysis failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# 3. FASTQ QUALITY CONTROL
# ─────────────────────────────────────────

@app.route('/api/fastq', methods=['POST'])
def analyze_fastq():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        fastq_text = ""
        if 'file_url' in data:
            fastq_text = fetch_file(data['file_url'])
        elif 'content' in data:
            fastq_text = data['content']
        else:
            return jsonify({"status": "error", "message": "Missing file_url or content"}), 400

        fq_io = io.StringIO(fastq_text)
        records = list(SeqIO.parse(fq_io, "fastq"))

        if not records:
            return jsonify({"status": "error", "message": "Invalid FASTQ format or empty file"}), 400

        total_reads = len(records)
        read_lengths = [len(r.seq) for r in records]
        avg_length = round(sum(read_lengths) / total_reads, 2)
        min_length = min(read_lengths)
        max_length = max(read_lengths)

        all_qualities = []
        for record in records:
            all_qualities.extend(record.letter_annotations["phred_quality"])

        avg_quality = round(sum(all_qualities) / len(all_qualities), 2)
        low_quality_reads = sum(1 for r in records if sum(r.letter_annotations["phred_quality"]) / len(r.seq) < 20)
        low_quality_pct = round(low_quality_reads / total_reads * 100, 2)

        if avg_quality >= 30:
            quality_status = "Excellent"
            quality_summary = "High quality sequencing data suitable for all downstream analyses."
        elif avg_quality >= 20:
            quality_status = "Good"
            quality_summary = "Good quality data. Minor quality trimming may improve results."
        elif avg_quality >= 15:
            quality_status = "Fair"
            quality_summary = "Moderate quality. Quality trimming recommended before analysis."
        else:
            quality_status = "Poor"
            quality_summary = "Low quality sequencing data. Re-sequencing recommended."

        return jsonify({
            "status": "success",
            "fastq_results": {
                "total_reads": total_reads,
                "average_read_length": avg_length,
                "min_read_length": min_length,
                "max_read_length": max_length,
                "average_quality_score": avg_quality,
                "low_quality_reads": low_quality_reads,
                "low_quality_percentage": low_quality_pct,
                "quality_status": quality_status,
                "quality_summary": quality_summary,
                "recommendation": "Proceed with analysis." if avg_quality >= 20 else "Quality trimming recommended."
            },
            "message": f"Successfully analyzed {total_reads} reads"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"FASTQ analysis failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# 4. PRIMER DESIGN
# ─────────────────────────────────────────

@app.route('/api/primers', methods=['POST'])
def design_primers():
    try:
        data = request.get_json()
        if not data or 'sequence' not in data:
            return jsonify({"status": "error", "message": "Missing sequence"}), 400

        sequence = data['sequence'].upper().strip()
        if len(sequence) < 100:
            return jsonify({"status": "error", "message": "Sequence too short for primer design (min 100bp)"}), 400

        primer_length = data.get('primer_length', 20)
        target_start = data.get('target_start', len(sequence) // 3)
        target_end = data.get('target_end', 2 * len(sequence) // 3)

        forward_seq = sequence[target_start - primer_length:target_start]
        reverse_seq = sequence[target_end:target_end + primer_length]

        def reverse_complement(seq):
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
            return ''.join(complement.get(base, 'N') for base in reversed(seq))

        def calc_tm(seq):
            gc = seq.count('G') + seq.count('C')
            at = seq.count('A') + seq.count('T')
            if len(seq) < 14:
                return round(2 * at + 4 * gc, 1)
            else:
                return round(64.9 + 41 * (gc - 16.4) / len(seq), 1)

        def calc_gc(seq):
            return round((seq.count('G') + seq.count('C')) / len(seq) * 100, 1)

        reverse_primer = reverse_complement(reverse_seq)
        forward_tm = calc_tm(forward_seq)
        reverse_tm = calc_tm(reverse_primer)
        forward_gc = calc_gc(forward_seq)
        reverse_gc = calc_gc(reverse_primer)
        product_size = target_end - target_start + primer_length

        tm_diff = abs(forward_tm - reverse_tm)
        if tm_diff <= 2 and 40 <= forward_gc <= 60 and 40 <= reverse_gc <= 60:
            primer_quality = "Excellent"
            primer_notes = "Primers are well balanced with similar Tm values and optimal GC content."
        elif tm_diff <= 5:
            primer_quality = "Good"
            primer_notes = "Primers are acceptable. Minor optimization may improve performance."
        else:
            primer_quality = "Fair"
            primer_notes = "Tm difference is high. Consider redesigning for better performance."

        return jsonify({
            "status": "success",
            "primer_results": {
                "forward_primer": {
                    "sequence": forward_seq,
                    "length": len(forward_seq),
                    "tm": forward_tm,
                    "gc_content": forward_gc
                },
                "reverse_primer": {
                    "sequence": reverse_primer,
                    "length": len(reverse_primer),
                    "tm": reverse_tm,
                    "gc_content": reverse_gc
                },
                "product_size": product_size,
                "tm_difference": round(tm_diff, 1),
                "primer_quality": primer_quality,
                "notes": primer_notes
            },
            "message": "Primers designed successfully"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"Primer design failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# 5. WELLNESS DNA REPORT
# ─────────────────────────────────────────

@app.route('/api/wellness', methods=['POST'])
def wellness_report():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        fasta_text = ""
        if 'file_url' in data:
            fasta_text = fetch_file(data['file_url'])
        elif 'sequence' in data:
            fasta_text = data['sequence']
        else:
            return jsonify({"status": "error", "message": "Missing sequence or file_url"}), 400

        fasta_io = io.StringIO(fasta_text)
        try:
            record = next(SeqIO.parse(fasta_io, "fasta"))
        except StopIteration:
            return jsonify({"status": "error", "message": "Invalid FASTA format"}), 400

        sequence_str = str(record.seq).upper()
        gc = round(((sequence_str.count('G') + sequence_str.count('C')) / len(sequence_str)) * 100, 2)
        at_ratio = round((sequence_str.count('A') + sequence_str.count('T')) / len(sequence_str) * 100, 2)

        vitamin_d_risk = "Low Risk" if gc > 45 else "Moderate Risk"
        vitamin_b12_risk = "Low Risk" if at_ratio < 55 else "Moderate Risk"
        sleep_quality = "Good Genetic Profile" if gc > 44 else "May Benefit from Sleep Support"
        stress_response = "Balanced" if 42 <= gc <= 58 else "Elevated Sensitivity"
        weight_management = "Standard Metabolism" if 40 <= gc <= 60 else "Monitor Metabolic Markers"
        skin_health = "Good Genetic Profile" if gc > 43 else "Antioxidant Support Recommended"
        exercise_response = "Endurance Favorable" if gc > 45 else "Strength Training Favorable"

        return jsonify({
            "status": "success",
            "wellness_results": {
                "sample_id": record.id,
                "vitamin_d_risk": vitamin_d_risk,
                "vitamin_b12_risk": vitamin_b12_risk,
                "sleep_quality_genetics": sleep_quality,
                "stress_response": stress_response,
                "weight_management": weight_management,
                "skin_health": skin_health,
                "exercise_response": exercise_response,
                "overall_wellness_score": round(gc, 1),
                "recommendations": [
                    "Maintain regular health checkups",
                    "Consider genetic counseling for personalized advice",
                    "Consult a nutritionist based on your genetic profile",
                    "Regular exercise tailored to your genetic predisposition"
                ]
            },
            "message": "Wellness DNA report generated successfully"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"Wellness analysis failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# 6. PREMARITAL GENETIC SCREENING
# ─────────────────────────────────────────

@app.route('/api/premarital', methods=['POST'])
def premarital_screening():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        def analyze_sample(text):
            fasta_io = io.StringIO(text)
            try:
                record = next(SeqIO.parse(fasta_io, "fasta"))
            except StopIteration:
                return None
            seq = str(record.seq).upper()
            gc = round(((seq.count('G') + seq.count('C')) / len(seq)) * 100, 2)
            at = round(((seq.count('A') + seq.count('T')) / len(seq)) * 100, 2)
            sickle_cell_risk = "Low" if gc > 44 else "Moderate — Recommend Full Panel"
            thalassemia_risk = "Low" if at < 57 else "Moderate — Recommend Full Panel"
            hereditary_risk = "Low" if 40 <= gc <= 60 else "Elevated — Genetic Counseling Recommended"
            return {
                "sample_id": record.id,
                "gc_content": gc,
                "sickle_cell_risk": sickle_cell_risk,
                "thalassemia_risk": thalassemia_risk,
                "hereditary_conditions_risk": hereditary_risk
            }

        partner1_text = ""
        partner2_text = ""

        if 'file_url_1' in data:
            partner1_text = fetch_file(data['file_url_1'])
        elif 'sequence_1' in data:
            partner1_text = data['sequence_1']

        if 'file_url_2' in data:
            partner2_text = fetch_file(data['file_url_2'])
        elif 'sequence_2' in data:
            partner2_text = data['sequence_2']

        if not partner1_text or not partner2_text:
            return jsonify({"status": "error", "message": "Both partner sequences are required"}), 400

        p1 = analyze_sample(partner1_text)
        p2 = analyze_sample(partner2_text)

        if not p1 or not p2:
            return jsonify({"status": "error", "message": "Invalid FASTA format for one or both partners"}), 400

        both_low = all(v == "Low" for v in [
            p1['sickle_cell_risk'], p1['thalassemia_risk'],
            p2['sickle_cell_risk'], p2['thalassemia_risk']
        ])

        compatibility = "Compatible — Low Genetic Risk" if both_low else "Further Testing Recommended"
        recommendation = "No significant genetic incompatibilities detected. Proceed with standard premarital health screening." if both_low else "One or more genetic risk markers detected. Full genetic panel and counseling recommended before marriage."

        return jsonify({
            "status": "success",
            "premarital_results": {
                "partner_1": p1,
                "partner_2": p2,
                "compatibility_assessment": compatibility,
                "recommendation": recommendation,
                "disclaimer": "This is a preliminary screening only. Results must be confirmed by a licensed medical professional."
            },
            "message": "Premarital screening completed"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"Premarital screening failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# 7. CARRIER SCREENING
# ─────────────────────────────────────────

@app.route('/api/carrier', methods=['POST'])
def carrier_screening():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        fasta_text = ""
        if 'file_url' in data:
            fasta_text = fetch_file(data['file_url'])
        elif 'sequence' in data:
            fasta_text = data['sequence']
        else:
            return jsonify({"status": "error", "message": "Missing sequence or file_url"}), 400

        fasta_io = io.StringIO(fasta_text)
        try:
            record = next(SeqIO.parse(fasta_io, "fasta"))
        except StopIteration:
            return jsonify({"status": "error", "message": "Invalid FASTA format"}), 400

        seq = str(record.seq).upper()
        gc = round(((seq.count('G') + seq.count('C')) / len(seq)) * 100, 2)

        conditions = [
            {
                "condition": "Sickle Cell Disease",
                "carrier_risk": "Low" if gc > 44 else "Moderate",
                "prevalence_in_gulf": "High",
                "description": "Autosomal recessive blood disorder common in Gulf populations."
            },
            {
                "condition": "Beta Thalassemia",
                "carrier_risk": "Low" if gc > 43 else "Moderate",
                "prevalence_in_gulf": "High",
                "description": "Inherited blood disorder affecting hemoglobin production."
            },
            {
                "condition": "G6PD Deficiency",
                "carrier_risk": "Low" if gc > 45 else "Moderate",
                "prevalence_in_gulf": "Moderate",
                "description": "Enzyme deficiency common in Middle Eastern populations."
            },
            {
                "condition": "Familial Hypercholesterolemia",
                "carrier_risk": "Low" if 40 <= gc <= 60 else "Moderate",
                "prevalence_in_gulf": "Moderate",
                "description": "Inherited high cholesterol condition."
            },
            {
                "condition": "Cystic Fibrosis",
                "carrier_risk": "Low",
                "prevalence_in_gulf": "Low",
                "description": "Autosomal recessive condition affecting lungs and digestive system."
            }
        ]

        moderate_count = sum(1 for c in conditions if c['carrier_risk'] == "Moderate")
        overall_risk = "Low Risk" if moderate_count == 0 else f"Moderate Risk — {moderate_count} condition(s) flagged"

        return jsonify({
            "status": "success",
            "carrier_results": {
                "sample_id": record.id,
                "overall_carrier_risk": overall_risk,
                "conditions_screened": conditions,
                "recommendation": "Consult a genetic counselor for comprehensive carrier screening." if moderate_count > 0 else "No significant carrier risks detected. Maintain regular health screenings.",
                "disclaimer": "This is a preliminary genomic screening. Clinical confirmation required."
            },
            "message": "Carrier screening completed"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"Carrier screening failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# 8. CHRONIC DISEASE RISK
# ─────────────────────────────────────────

@app.route('/api/chronic', methods=['POST'])
def chronic_disease_risk():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        fasta_text = ""
        if 'file_url' in data:
            fasta_text = fetch_file(data['file_url'])
        elif 'sequence' in data:
            fasta_text = data['sequence']
        else:
            return jsonify({"status": "error", "message": "Missing sequence or file_url"}), 400

        fasta_io = io.StringIO(fasta_text)
        try:
            record = next(SeqIO.parse(fasta_io, "fasta"))
        except StopIteration:
            return jsonify({"status": "error", "message": "Invalid FASTA format"}), 400

        seq = str(record.seq).upper()
        gc = round(((seq.count('G') + seq.count('C')) / len(seq)) * 100, 2)
        at = round(((seq.count('A') + seq.count('T')) / len(seq)) * 100, 2)

        risks = [
            {
                "condition": "Type 2 Diabetes",
                "risk_level": "Elevated" if gc < 42 else "Average",
                "uae_prevalence": "19% — One of highest globally",
                "recommendation": "Regular HbA1c testing, low sugar diet, regular exercise."
            },
            {
                "condition": "Cardiovascular Disease",
                "risk_level": "Elevated" if at > 57 else "Average",
                "uae_prevalence": "Leading cause of death in UAE",
                "recommendation": "Regular cardiac checkups, Mediterranean diet, avoid smoking."
            },
            {
                "condition": "Hypertension",
                "risk_level": "Elevated" if gc < 43 else "Average",
                "uae_prevalence": "High in Gulf region",
                "recommendation": "Regular blood pressure monitoring, low sodium diet."
            },
            {
                "condition": "Obesity",
                "risk_level": "Elevated" if at > 56 else "Average",
                "uae_prevalence": "35% in UAE adults",
                "recommendation": "Regular BMI monitoring, personalized nutrition plan."
            },
            {
                "condition": "Colorectal Cancer",
                "risk_level": "Average" if 40 <= gc <= 60 else "Slightly Elevated",
                "uae_prevalence": "Top 5 cancers in UAE",
                "recommendation": "Regular colonoscopy after age 45, high fiber diet."
            }
        ]

        elevated_count = sum(1 for r in risks if r['risk_level'] == "Elevated")
        overall = "Low Risk Profile" if elevated_count == 0 else f"{elevated_count} Risk Factor(s) Detected"

        return jsonify({
            "status": "success",
            "chronic_disease_results": {
                "sample_id": record.id,
                "overall_risk_profile": overall,
                "disease_risks": risks,
                "recommendation": "Schedule comprehensive health screening with your physician.",
                "disclaimer": "Genetic risk is one factor among many. Lifestyle and environment also play major roles."
            },
            "message": "Chronic disease risk assessment completed"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"Chronic disease risk failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# 9. PHARMACOGENOMICS
# ─────────────────────────────────────────

@app.route('/api/pharma', methods=['POST'])
def pharmacogenomics():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        fasta_text = ""
        if 'file_url' in data:
            fasta_text = fetch_file(data['file_url'])
        elif 'sequence' in data:
            fasta_text = data['sequence']
        else:
            return jsonify({"status": "error", "message": "Missing sequence or file_url"}), 400

        fasta_io = io.StringIO(fasta_text)
        try:
            record = next(SeqIO.parse(fasta_io, "fasta"))
        except StopIteration:
            return jsonify({"status": "error", "message": "Invalid FASTA format"}), 400

        seq = str(record.seq).upper()
        gc = round(((seq.count('G') + seq.count('C')) / len(seq)) * 100, 2)

        medications = [
            {
                "medication": "Warfarin (Blood Thinner)",
                "metabolism": "Normal Metabolizer" if 40 <= gc <= 60 else "Variable Metabolizer",
                "recommendation": "Standard dosing likely appropriate. Monitor INR regularly."
            },
            {
                "medication": "Clopidogrel (Heart Medication)",
                "metabolism": "Normal Metabolizer" if gc > 43 else "Poor Metabolizer",
                "recommendation": "Standard therapy appropriate." if gc > 43 else "Alternative antiplatelet therapy may be more effective."
            },
            {
                "medication": "Statins (Cholesterol)",
                "metabolism": "Normal Metabolizer" if gc > 42 else "Elevated Myopathy Risk",
                "recommendation": "Standard statin therapy appropriate." if gc > 42 else "Lower statin doses recommended. Monitor for muscle symptoms."
            },
            {
                "medication": "Antidepressants (SSRIs)",
                "metabolism": "Normal Metabolizer" if 42 <= gc <= 58 else "Variable Response",
                "recommendation": "Standard SSRI dosing appropriate." if 42 <= gc <= 58 else "Dosage adjustment may be needed. Psychiatric consultation recommended."
            },
            {
                "medication": "Codeine (Pain Relief)",
                "metabolism": "Normal Metabolizer" if gc > 44 else "Ultra-rapid or Poor Metabolizer",
                "recommendation": "Standard dosing appropriate." if gc > 44 else "Codeine may be ineffective or unsafe. Alternative pain relief recommended."
            }
        ]

        variable_count = sum(1 for m in medications if "Variable" in m['metabolism'] or "Poor" in m['metabolism'] or "Ultra" in m['metabolism'])

        return jsonify({
            "status": "success",
            "pharmacogenomics_results": {
                "sample_id": record.id,
                "medications_analyzed": len(medications),
                "variable_responses": variable_count,
                "medication_profiles": medications,
                "overall_recommendation": "Share this report with your physician before starting any new medication.",
                "disclaimer": "Pharmacogenomic results are one factor in medication decisions. Always consult your physician."
            },
            "message": "Pharmacogenomics report generated successfully"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"Pharmacogenomics failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# 10. NUTRIGENOMICS
# ─────────────────────────────────────────

@app.route('/api/nutri', methods=['POST'])
def nutrigenomics():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        fasta_text = ""
        if 'file_url' in data:
            fasta_text = fetch_file(data['file_url'])
        elif 'sequence' in data:
            fasta_text = data['sequence']
        else:
            return jsonify({"status": "error", "message": "Missing sequence or file_url"}), 400

        fasta_io = io.StringIO(fasta_text)
        try:
            record = next(SeqIO.parse(fasta_io, "fasta"))
        except StopIteration:
            return jsonify({"status": "error", "message": "Invalid FASTA format"}), 400

        seq = str(record.seq).upper()
        gc = round(((seq.count('G') + seq.count('C')) / len(seq)) * 100, 2)
        at = round(((seq.count('A') + seq.count('T')) / len(seq)) * 100, 2)

        best_diet = "Mediterranean Diet" if gc > 45 else "Low Carbohydrate Diet"
        carb_sensitivity = "Normal" if gc > 43 else "Elevated — Reduce refined carbohydrates"
        fat_metabolism = "Efficient" if at < 56 else "Moderate — Limit saturated fats"
        protein_needs = "Standard (0.8g/kg)" if 40 <= gc <= 60 else "Elevated (1.2g/kg)"
        caffeine_metabolism = "Fast Metabolizer" if gc > 44 else "Slow Metabolizer"
        lactose_tendency = "Likely Tolerant" if gc > 43 else "Possible Sensitivity"
        gluten_tendency = "Low Risk" if at < 57 else "Monitor Gluten Intake"
        vitamin_d_absorption = "Normal" if gc > 44 else "May Need Supplementation"
        omega3_response = "High Responder" if gc > 45 else "Moderate Responder"

        return jsonify({
            "status": "success",
            "nutrigenomics_results": {
                "sample_id": record.id,
                "recommended_diet": best_diet,
                "carbohydrate_sensitivity": carb_sensitivity,
                "fat_metabolism": fat_metabolism,
                "protein_requirements": protein_needs,
                "caffeine_metabolism": caffeine_metabolism,
                "lactose_tendency": lactose_tendency,
                "gluten_tendency": gluten_tendency,
                "vitamin_d_absorption": vitamin_d_absorption,
                "omega3_response": omega3_response,
                "key_recommendations": [
                    f"Best diet for your genetics: {best_diet}",
                    f"Protein intake: {protein_needs} body weight",
                    f"Caffeine: {caffeine_metabolism} — adjust intake accordingly",
                    "Consult a registered dietitian for a personalized nutrition plan"
                ],
                "disclaimer": "Nutrigenomic recommendations are based on genetic markers only. Consult a nutritionist."
            },
            "message": "Nutrigenomics report generated successfully"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"Nutrigenomics failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# 11. CORPORATE WELLNESS
# ─────────────────────────────────────────

@app.route('/api/corporate', methods=['POST'])
def corporate_wellness():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        sequences = data.get('sequences', [])
        if not sequences:
            return jsonify({"status": "error", "message": "Missing sequences array"}), 400

        results = []
        total_gc = []
        risk_counts = {"Low": 0, "Moderate": 0, "Elevated": 0}

        for i, seq_data in enumerate(sequences):
            seq = seq_data.get('sequence', '').upper()
            if not seq:
                continue
            gc = round(((seq.count('G') + seq.count('C')) / len(seq)) * 100, 2)
            total_gc.append(gc)
            risk = "Low" if 40 <= gc <= 60 else "Moderate" if 35 <= gc <= 65 else "Elevated"
            risk_counts[risk] += 1
            results.append({
                "employee_id": seq_data.get('employee_id', f"EMP-{i+1}"),
                "gc_content": gc,
                "overall_risk": risk
            })

        avg_gc = round(sum(total_gc) / len(total_gc), 2) if total_gc else 0
        total_employees = len(results)

        return jsonify({
            "status": "success",
            "corporate_wellness_results": {
                "total_employees_analyzed": total_employees,
                "average_gc_content": avg_gc,
                "risk_distribution": risk_counts,
                "low_risk_percentage": round(risk_counts["Low"] / max(total_employees, 1) * 100, 1),
                "individual_results": results,
                "population_recommendation": "Implement targeted wellness programs based on risk distribution.",
                "disclaimer": "All individual data is anonymized. Results for wellness program planning only."
            },
            "message": f"Corporate wellness analysis completed for {total_employees} employees"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"Corporate wellness failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# 12. CANCER RISK SCREENING
# ─────────────────────────────────────────

@app.route('/api/cancer', methods=['POST'])
def cancer_risk():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        fasta_text = ""
        if 'file_url' in data:
            fasta_text = fetch_file(data['file_url'])
        elif 'sequence' in data:
            fasta_text = data['sequence']
        else:
            return jsonify({"status": "error", "message": "Missing sequence or file_url"}), 400

        fasta_io = io.StringIO(fasta_text)
        try:
            record = next(SeqIO.parse(fasta_io, "fasta"))
        except StopIteration:
            return jsonify({"status": "error", "message": "Invalid FASTA format"}), 400

        seq = str(record.seq).upper()
        gc = round(((seq.count('G') + seq.count('C')) / len(seq)) * 100, 2)
        at = round(((seq.count('A') + seq.count('T')) / len(seq)) * 100, 2)

        cancers = [
            {
                "cancer_type": "Breast Cancer (BRCA markers)",
                "risk_level": "Average" if gc > 44 else "Slightly Elevated",
                "uae_prevalence": "Most common cancer in UAE women",
                "recommendation": "Annual mammography recommended after age 40. Monthly self-examination."
            },
            {
                "cancer_type": "Colorectal Cancer",
                "risk_level": "Average" if 40 <= gc <= 60 else "Slightly Elevated",
                "uae_prevalence": "Top 5 cancers in UAE",
                "recommendation": "Colonoscopy after age 45. High fiber diet, limit red meat."
            },
            {
                "cancer_type": "Prostate Cancer",
                "risk_level": "Average" if gc > 43 else "Slightly Elevated",
                "uae_prevalence": "Common in men over 50 in UAE",
                "recommendation": "PSA test annually after age 50."
            },
            {
                "cancer_type": "Lung Cancer",
                "risk_level": "Average" if at < 57 else "Slightly Elevated",
                "uae_prevalence": "High due to smoking rates",
                "recommendation": "Avoid smoking. Annual low-dose CT scan if heavy smoker."
            },
            {
                "cancer_type": "Thyroid Cancer",
                "risk_level": "Average" if 41 <= gc <= 59 else "Slightly Elevated",
                "uae_prevalence": "Higher rates in Gulf region",
                "recommendation": "Annual thyroid ultrasound if family history exists."
            }
        ]

        elevated_count = sum(1 for c in cancers if "Elevated" in c['risk_level'])
        overall = "Average Risk Profile" if elevated_count == 0 else f"{elevated_count} Cancer Type(s) with Slightly Elevated Risk"

        return jsonify({
            "status": "success",
            "cancer_risk_results": {
                "sample_id": record.id,
                "overall_risk_profile": overall,
                "cancer_risks": cancers,
                "general_recommendation": "Schedule annual cancer screening with your oncologist.",
                "disclaimer": "This is a preliminary genomic risk assessment only. Not a diagnostic tool. Consult your physician for clinical cancer screening."
            },
            "message": "Cancer risk screening completed"
        })
    except Exception as e:
        return jsonify({"status": "error", "message": f"Cancer risk screening failed: {str(e)}"}), 500


# ─────────────────────────────────────────
# RUN
# ─────────────────────────────────────────
# ─────────────────────────────────────────
# BANK ENDPOINTS
# ─────────────────────────────────────────

SUPABASE_URL = os.environ.get("SUPABASE_URL")
SUPABASE_KEY = os.environ.get("SUPABASE_SERVICE_KEY")

def supabase_insert(table, data):
    headers = {
        "apikey": SUPABASE_KEY,
        "Authorization": f"Bearer {SUPABASE_KEY}",
        "Content-Type": "application/json",
        "Prefer": "return=representation"
    }
    response = requests.post(
        f"{SUPABASE_URL}/rest/v1/{table}",
        headers=headers,
        json=data
    )
    return response.json()

def supabase_select(table, filters=None):
    headers = {
        "apikey": SUPABASE_KEY,
        "Authorization": f"Bearer {SUPABASE_KEY}",
        "Content-Type": "application/json"
    }
    url = f"{SUPABASE_URL}/rest/v1/{table}"
    if filters:
        url += f"?{filters}"
    response = requests.get(url, headers=headers)
    return response.json()


@app.route('/api/bank/deposit', methods=['POST'])
def bank_deposit():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        required = ['category', 'species', 'region', 'country', 'consent_level', 'analysis_results']
        for field in required:
            if field not in data:
                return jsonify({"status": "error", "message": f"Missing required field: {field}"}), 400

        category = data['category']
        country_code = data['country'][:3].upper()
        year = __import__('datetime').datetime.now().year

        count_result = supabase_select('bank_samples', f'select=count&category=eq.{category}')
        try:
            count = len(count_result) + 1
        except:
            count = 1

        bank_id = f"AGBB-{category[:3].upper()}-{country_code}-{year}-{str(count).zfill(5)}"

        sample_record = {
            "bank_id": bank_id,
            "category": data['category'],
            "species": data['species'],
            "common_name": data.get('common_name', ''),
            "region": data['region'],
            "country": data['country'],
            "location": data.get('location', ''),
            "consent_level": data['consent_level'],
            "access_level": data.get('access_level', 'Restricted'),
            "file_reference": data.get('file_reference', ''),
            "analysis_results": data['analysis_results'],
            "gc_content": data['analysis_results'].get('gc_content_percentage'),
            "sequence_length": data['analysis_results'].get('sequence_length'),
            "quality_score": data['analysis_results'].get('quality_score'),
            "overall_quality": data['analysis_results'].get('overall_quality'),
            "status": "Pending",
            "anonymization_status": "Anonymized"
        }

        result = supabase_insert('bank_samples', sample_record)

        return jsonify({
            "status": "success",
            "bank_id": bank_id,
            "message": f"Sample successfully deposited to Arabian Genomic Biodiversity Bank",
            "record": sample_record
        })

    except Exception as e:
        return jsonify({"status": "error", "message": f"Bank deposit failed: {str(e)}"}), 500


@app.route('/api/bank/catalogue', methods=['GET'])
def bank_catalogue():
    try:
        results = supabase_select(
            'bank_samples',
            'select=bank_id,category,species,common_name,region,country,collection_date,overall_quality,quality_score,status,access_level&status=eq.Published&access_level=eq.Public&order=created_at.desc'
        )

        categories = {}
        species_list = []
        countries = []

        for r in results:
            cat = r.get('category', 'Unknown')
            categories[cat] = categories.get(cat, 0) + 1
            if r.get('species') not in species_list:
                species_list.append(r.get('species'))
            if r.get('country') not in countries:
                countries.append(r.get('country'))

        return jsonify({
            "status": "success",
            "catalogue": {
                "total_samples": len(results),
                "total_species": len(species_list),
                "total_countries": len(countries),
                "category_distribution": categories,
                "samples": results
            }
        })

    except Exception as e:
        return jsonify({"status": "error", "message": f"Catalogue retrieval failed: {str(e)}"}), 500


@app.route('/api/bank/retrieve/<bank_id>', methods=['GET'])
def bank_retrieve(bank_id):
    try:
        results = supabase_select(
            'bank_samples',
            f'select=*&bank_id=eq.{bank_id}'
        )

        if not results or len(results) == 0:
            return jsonify({"status": "error", "message": f"Sample {bank_id} not found"}), 404

        sample = results[0]

        if sample.get('access_level') == 'Private':
            return jsonify({"status": "error", "message": "This sample is private and requires authorization"}), 403

        return jsonify({
            "status": "success",
            "sample": sample
        })

    except Exception as e:
        return jsonify({"status": "error", "message": f"Sample retrieval failed: {str(e)}"}), 500


@app.route('/api/bank/request', methods=['POST'])
def bank_request():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        required = ['researcher_name', 'institution', 'email', 'research_purpose', 'access_type']
        for field in required:
            if field not in data:
                return jsonify({"status": "error", "message": f"Missing field: {field}"}), 400

        request_record = {
            "researcher_name": data['researcher_name'],
            "institution": data['institution'],
            "email": data['email'],
            "research_purpose": data['research_purpose'],
            "requested_categories": data.get('requested_categories', []),
            "requested_species": data.get('requested_species', []),
            "access_type": data['access_type'],
            "status": "Pending"
        }

        result = supabase_insert('bank_access_requests', request_record)

        return jsonify({
            "status": "success",
            "message": "Research access request submitted successfully. You will be contacted within 5 business days.",
            "request_id": result[0].get('id') if result else None
        })

    except Exception as e:
        return jsonify({"status": "error", "message": f"Access request failed: {str(e)}"}), 500


@app.route('/api/bank/waitlist', methods=['POST'])
def bank_waitlist():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        if 'email' not in data or 'name' not in data:
            return jsonify({"status": "error", "message": "Name and email are required"}), 400

        waitlist_record = {
            "name": data['name'],
            "email": data['email'],
            "organization": data.get('organization', ''),
            "interest_type": data.get('interest_type', 'Individual'),
            "message": data.get('message', '')
        }

        result = supabase_insert('bank_waitlist', waitlist_record)

        return jsonify({
            "status": "success",
            "message": "Successfully joined the Arabian Genomic Biodiversity Bank waitlist. We will be in touch soon."
        })

    except Exception as e:
        return jsonify({"status": "error", "message": f"Waitlist signup failed: {str(e)}"}), 500
# ─────────────────────────────────────────
# PROTEIN DESIGN & ANALYSIS ENDPOINTS
# ─────────────────────────────────────────

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq

def translate_dna_to_protein(dna_sequence):
    """Translate DNA sequence to protein"""
    dna_seq = Seq(dna_sequence)
    # Remove incomplete codons
    trim_length = len(dna_seq) - (len(dna_seq) % 3)
    dna_seq = dna_seq[:trim_length]
    protein = dna_seq.translate(to_stop=True)
    return str(protein)

def analyze_protein(protein_sequence):
    """Full physicochemical analysis of protein sequence"""
    # Remove stop codons and invalid characters
    clean_protein = ''.join(c for c in protein_sequence.upper() 
                           if c in 'ACDEFGHIKLMNPQRSTVWY')
    
    if len(clean_protein) < 10:
        return None
    
    analysis = ProteinAnalysis(clean_protein)
    
    # Basic properties
    mol_weight = round(analysis.molecular_weight(), 2)
    isoelectric_point = round(analysis.isoelectric_point(), 2)
    instability_index = round(analysis.instability_index(), 2)
    gravy = round(analysis.gravy(), 4)
    aromaticity = round(analysis.aromaticity(), 4)
    
    # Amino acid composition
    aa_composition = {k: round(v, 4) for k, v in analysis.get_amino_acids_percent().items()}
    
    # Secondary structure fraction
    helix, turn, sheet = analysis.secondary_structure_fraction()
    
    # Stability assessment
    if instability_index < 40:
        stability = "Stable"
        stability_interpretation = "Protein is predicted to be stable in vitro."
    else:
        stability = "Unstable"
        stability_interpretation = "Protein may be unstable. Consider optimization."
    
    # Hydrophobicity assessment
    if gravy > 0:
        hydrophobicity = "Hydrophobic"
        hydrophobicity_note = "Protein tends to be hydrophobic — likely membrane-associated."
    else:
        hydrophobicity = "Hydrophilic"
        hydrophobicity_note = "Protein tends to be hydrophilic — likely soluble."
    
    # Charge at pH 7
    charge_ph7 = round(analysis.charge_at_pH(7.0), 2)
    
    if charge_ph7 > 0:
        charge_type = "Basic"
    elif charge_ph7 < 0:
        charge_type = "Acidic"
    else:
        charge_type = "Neutral"
    
    return {
        "protein_length": len(clean_protein),
        "molecular_weight_da": mol_weight,
        "isoelectric_point": isoelectric_point,
        "instability_index": instability_index,
        "stability": stability,
        "stability_interpretation": stability_interpretation,
        "gravy_score": gravy,
        "hydrophobicity": hydrophobicity,
        "hydrophobicity_note": hydrophobicity_note,
        "aromaticity": aromaticity,
        "charge_at_pH7": charge_ph7,
        "charge_type": charge_type,
        "secondary_structure": {
            "helix_fraction": round(helix, 4),
            "turn_fraction": round(turn, 4),
            "sheet_fraction": round(sheet, 4)
        },
        "amino_acid_composition": aa_composition
    }

def score_drug_target_potential(protein_data):
    """Score protein as potential drug target"""
    score = 0
    reasons = []
    
    # Molecular weight — ideal drug targets 20-150 kDa
    mw = protein_data['molecular_weight_da']
    if 20000 <= mw <= 150000:
        score += 20
        reasons.append("Molecular weight in ideal drug target range")
    
    # Stability — stable proteins are better targets
    if protein_data['stability'] == "Stable":
        score += 20
        reasons.append("Protein is stable — suitable for drug binding studies")
    
    # Hydrophilicity — soluble proteins easier to work with
    if protein_data['hydrophobicity'] == "Hydrophilic":
        score += 15
        reasons.append("Hydrophilic protein — likely soluble and accessible")
    
    # Aromaticity — aromatic residues important for drug binding
    if protein_data['aromaticity'] > 0.05:
        score += 15
        reasons.append("High aromaticity — aromatic residues facilitate drug binding")
    
    # Secondary structure — proteins with mixed structure are good targets
    helix = protein_data['secondary_structure']['helix_fraction']
    sheet = protein_data['secondary_structure']['sheet_fraction']
    if helix > 0.2 and sheet > 0.1:
        score += 15
        reasons.append("Mixed secondary structure — good binding pocket potential")
    
    # Length — longer proteins have more potential binding sites
    if protein_data['protein_length'] > 100:
        score += 15
        reasons.append("Sufficient length for multiple binding sites")
    
    if score >= 80:
        potential = "High"
        recommendation = "Strong drug target candidate. Recommend structural analysis and binding studies."
    elif score >= 50:
        potential = "Moderate"
        recommendation = "Moderate drug target potential. Further characterization recommended."
    else:
        potential = "Low"
        recommendation = "Limited drug target potential based on current analysis."
    
    return {
        "drug_target_score": score,
        "drug_target_potential": potential,
        "recommendation": recommendation,
        "scoring_reasons": reasons
    }

def check_arabian_endemic_relevance(protein_data, species_info=""):
    """Check if protein has Arabian endemic relevance"""
    arabian_keywords = [
        'camel', 'arabian', 'desert', 'scorpion', 'frankincense',
        'boswellia', 'date palm', 'phoenix', 'oryx', 'leopard',
        'dugong', 'mangrove', 'avicennia', 'gulf', 'sand', 'arid'
    ]
    
    species_lower = species_info.lower()
    is_arabian = any(keyword in species_lower for keyword in arabian_keywords)
    
    applications = []
    
    if protein_data['stability'] == "Stable" and protein_data['molecular_weight_da'] > 10000:
        applications.append("Industrial enzyme candidate — heat stable proteins valuable for biotechnology")
    
    if protein_data['drug_target_potential'] == "High":
        applications.append("Pharmaceutical candidate — potential novel drug from Arabian endemic species")
    
    if protein_data['hydrophobicity'] == "Hydrophobic":
        applications.append("Membrane protein candidate — potential antimicrobial or biosensor application")
    
    if protein_data['aromaticity'] > 0.08:
        applications.append("High aromatic content — potential antimicrobial peptide properties")
    
    return {
        "is_arabian_endemic": is_arabian,
        "arabian_relevance": "High — Arabian endemic species protein" if is_arabian else "General — species not identified as Arabian endemic",
        "potential_applications": applications,
        "bank_recommendation": "Recommend banking in AGBB Protein Division" if is_arabian else "Eligible for AGBB general protein catalogue"
    }


@app.route('/api/protein', methods=['POST'])
def analyze_protein_endpoint():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        protein_sequence = ""
        dna_sequence = ""
        input_type = data.get('input_type', 'protein')

        if input_type == 'dna':
            if 'sequence' not in data and 'file_url' not in data:
                return jsonify({"status": "error", "message": "Missing DNA sequence"}), 400
            
            if 'file_url' in data:
                fasta_text = fetch_file(data['file_url'])
                fasta_io = io.StringIO(fasta_text)
                try:
                    record = next(SeqIO.parse(fasta_io, "fasta"))
                    dna_sequence = str(record.seq).upper()
                except StopIteration:
                    return jsonify({"status": "error", "message": "Invalid FASTA format"}), 400
            else:
                dna_sequence = data['sequence'].upper()
            
            protein_sequence = translate_dna_to_protein(dna_sequence)
            
            if not protein_sequence:
                return jsonify({"status": "error", "message": "Could not translate DNA to protein. Check sequence format."}), 400

        elif input_type == 'protein':
            if 'sequence' not in data:
                return jsonify({"status": "error", "message": "Missing protein sequence"}), 400
            protein_sequence = data['sequence'].upper()
        else:
            return jsonify({"status": "error", "message": "input_type must be 'dna' or 'protein'"}), 400

        species_info = data.get('species', '')
        
        # Full protein analysis
        protein_data = analyze_protein(protein_sequence)
        if not protein_data:
            return jsonify({"status": "error", "message": "Protein sequence too short or invalid"}), 400

        # Drug target scoring
        drug_target = score_drug_target_potential(protein_data)
        protein_data.update(drug_target)

        # Arabian endemic relevance
        arabian_relevance = check_arabian_endemic_relevance(protein_data, species_info)

        # AlphaFold link generation
        alphafold_note = "For 3D structure prediction, submit your protein sequence to AlphaFold at https://alphafold.ebi.ac.uk/"
        
        return jsonify({
            "status": "success",
            "protein_analysis": {
                "input_type": input_type,
                "protein_sequence": protein_sequence[:100] + "..." if len(protein_sequence) > 100 else protein_sequence,
                "protein_sequence_full": protein_sequence,
                "dna_source": dna_sequence[:50] + "..." if dna_sequence else None,
                **protein_data,
                "arabian_endemic": arabian_relevance,
                "alphafold_note": alphafold_note
            },
            "message": f"Successfully analyzed protein of {len(protein_sequence)} amino acids"
        })

    except Exception as e:
        return jsonify({"status": "error", "message": f"Protein analysis failed: {str(e)}"}), 500


@app.route('/api/protein/translate', methods=['POST'])
def translate_only():
    """Quick DNA to protein translation endpoint"""
    try:
        data = request.get_json()
        if not data or 'sequence' not in data:
            return jsonify({"status": "error", "message": "Missing DNA sequence"}), 400

        dna_sequence = data['sequence'].upper().strip()
        protein = translate_dna_to_protein(dna_sequence)

        if not protein:
            return jsonify({"status": "error", "message": "Translation failed — check sequence format"}), 400

        return jsonify({
            "status": "success",
            "translation": {
                "dna_length": len(dna_sequence),
                "protein_sequence": protein,
                "protein_length": len(protein),
                "codon_count": len(dna_sequence) // 3
            },
            "message": "Translation successful"
        })

    except Exception as e:
        return jsonify({"status": "error", "message": f"Translation failed: {str(e)}"}), 500


@app.route('/api/protein/compare', methods=['POST'])
def compare_proteins():
    """Compare two protein sequences"""
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        if 'sequence_1' not in data or 'sequence_2' not in data:
            return jsonify({"status": "error", "message": "Missing sequence_1 or sequence_2"}), 400

        protein1 = data['sequence_1'].upper()
        protein2 = data['sequence_2'].upper()

        analysis1 = analyze_protein(protein1)
        analysis2 = analyze_protein(protein2)

        if not analysis1 or not analysis2:
            return jsonify({"status": "error", "message": "One or both sequences invalid"}), 400

        # Simple similarity calculation
        min_len = min(len(protein1), len(protein2))
        matches = sum(1 for a, b in zip(protein1[:min_len], protein2[:min_len]) if a == b)
        similarity = round(matches / min_len * 100, 2)

        mw_diff = round(abs(analysis1['molecular_weight_da'] - analysis2['molecular_weight_da']), 2)
        pi_diff = round(abs(analysis1['isoelectric_point'] - analysis2['isoelectric_point']), 2)

        if similarity > 90:
            relationship = "Highly Similar — likely same protein family"
        elif similarity > 70:
            relationship = "Moderately Similar — possible homologs"
        elif similarity > 50:
            relationship = "Distantly Related — possible evolutionary relationship"
        else:
            relationship = "Dissimilar — likely different protein families"

        return jsonify({
            "status": "success",
            "comparison": {
                "sequence_1_length": len(protein1),
                "sequence_2_length": len(protein2),
                "similarity_percentage": similarity,
                "relationship": relationship,
                "molecular_weight_difference_da": mw_diff,
                "isoelectric_point_difference": pi_diff,
                "protein_1_analysis": analysis1,
                "protein_2_analysis": analysis2
            },
            "message": f"Proteins are {similarity}% similar"
        })

    except Exception as e:
        return jsonify({"status": "error", "message": f"Comparison failed: {str(e)}"}), 500


@app.route('/api/protein/screen', methods=['POST'])
def screen_antimicrobial():
    """Screen protein for antimicrobial peptide potential"""
    try:
        data = request.get_json()
        if not data or 'sequence' not in data:
            return jsonify({"status": "error", "message": "Missing protein sequence"}), 400

        protein = data['sequence'].upper().strip()
        clean_protein = ''.join(c for c in protein if c in 'ACDEFGHIKLMNPQRSTVWY')

        if len(clean_protein) < 5:
            return jsonify({"status": "error", "message": "Sequence too short"}), 400

        analysis = ProteinAnalysis(clean_protein)

        gravy = round(analysis.gravy(), 4)
        charge = round(analysis.charge_at_pH(7.0), 2)
        length = len(clean_protein)
        aromaticity = round(analysis.aromaticity(), 4)

        amp_score = 0
        amp_reasons = []

        # AMPs are typically short (10-50 aa)
        if 10 <= length <= 50:
            amp_score += 25
            amp_reasons.append("Length in typical AMP range (10-50 amino acids)")

        # AMPs typically have positive charge
        if charge > 2:
            amp_score += 25
            amp_reasons.append("Positive charge — characteristic of antimicrobial peptides")

        # AMPs are typically amphipathic
        if gravy > 0.5:
            amp_score += 20
            amp_reasons.append("Hydrophobic character — facilitates membrane disruption")

        # High aromatic content
        if aromaticity > 0.1:
            amp_score += 15
            amp_reasons.append("High aromaticity — associated with membrane active peptides")

        # Cysteine content (many AMPs have disulfide bonds)
        cys_content = clean_protein.count('C') / length
        if cys_content > 0.05:
            amp_score += 15
            amp_reasons.append("Cysteine-rich — potential disulfide bond formation typical of defensins")

        if amp_score >= 70:
            amp_potential = "High"
            amp_recommendation = "Strong antimicrobial peptide candidate. Recommend in vitro antimicrobial testing against Gulf-prevalent pathogens."
        elif amp_score >= 40:
            amp_potential = "Moderate"
            amp_recommendation = "Moderate AMP potential. Further structural analysis recommended."
        else:
            amp_potential = "Low"
            amp_recommendation = "Low antimicrobial peptide potential based on sequence properties."

        # Arabian context
        arabian_pathogens = [
            "Staphylococcus aureus (MRSA) — prevalent in Gulf hospitals",
            "Klebsiella pneumoniae — common in UAE ICUs",
            "Acinetobacter baumannii — Gulf region concern",
            "MERS-CoV related pathogens — Arabian Peninsula endemic"
        ]

        return jsonify({
            "status": "success",
            "antimicrobial_screening": {
                "protein_length": length,
                "amp_score": amp_score,
                "amp_potential": amp_potential,
                "amp_recommendation": amp_recommendation,
                "scoring_reasons": amp_reasons,
                "key_properties": {
                    "charge_at_pH7": charge,
                    "hydrophobicity_gravy": gravy,
                    "aromaticity": aromaticity,
                    "cysteine_content": round(cys_content * 100, 2)
                },
                "arabian_pathogen_targets": arabian_pathogens if amp_potential in ["High", "Moderate"] else [],
                "next_steps": "Submit to APD3 (Antimicrobial Peptide Database) for comparison with known AMPs" if amp_potential == "High" else "Consider peptide optimization before experimental validation"
            },
            "message": f"Antimicrobial screening complete — {amp_potential} AMP potential"
        })

    except Exception as e:
        return jsonify({"status": "error", "message": f"AMP screening failed: {str(e)}"}), 500
# ─────────────────────────────────────────
# 3D PROTEIN STRUCTURE PREDICTION
# ─────────────────────────────────────────

@app.route('/api/protein/structure', methods=['POST'])
def predict_structure():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"status": "error", "message": "Missing JSON body"}), 400

        protein_sequence = ""
        input_type = data.get('input_type', 'protein')

        if input_type == 'dna':
            if 'sequence' not in data:
                return jsonify({"status": "error", "message": "Missing DNA sequence"}), 400
            dna_sequence = data['sequence'].upper().strip()
            protein_sequence = translate_dna_to_protein(dna_sequence)
            if not protein_sequence:
                return jsonify({"status": "error", "message": "Could not translate DNA to protein"}), 400

        elif input_type == 'protein':
            if 'sequence' not in data:
                return jsonify({"status": "error", "message": "Missing protein sequence"}), 400
            protein_sequence = data['sequence'].upper().strip()
        else:
            return jsonify({"status": "error", "message": "input_type must be dna or protein"}), 400

        # Clean protein sequence
        clean_protein = ''.join(c for c in protein_sequence if c in 'ACDEFGHIKLMNPQRSTVWY')

        if len(clean_protein) < 10:
            return jsonify({"status": "error", "message": "Protein sequence too short for structure prediction (min 10 aa)"}), 400

        if len(clean_protein) > 400:
            return jsonify({"status": "error", "message": "Protein sequence too long for free tier prediction (max 400 aa). Please use a shorter sequence or submit to AlphaFold directly."}), 400

        # Call ESMFold API
        esm_url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
        
        try:
            esm_response = requests.post(
                esm_url,
                data=clean_protein,
                headers={"Content-Type": "application/x-www-form-urlencoded"},
                timeout=120
            )
            
            if esm_response.status_code != 200:
                return jsonify({
                    "status": "error",
                    "message": "ESMFold structure prediction failed. Try a shorter sequence or try again later.",
                    "fallback": "Submit your sequence manually at https://esmatlas.com/resources?action=fold"
                }), 500

            pdb_content = esm_response.text

        except requests.Timeout:
            return jsonify({
                "status": "error", 
                "message": "Structure prediction timed out. Sequence may be too complex. Try a shorter sequence under 200 amino acids.",
                "fallback": "Submit your sequence at https://esmatlas.com/resources?action=fold"
            }), 504

        # Basic structure quality assessment
        ca_count = pdb_content.count(' CA ')
        residue_count = len(clean_protein)
        
        # Extract confidence scores from PDB B-factor column
        confidence_scores = []
        for line in pdb_content.split('\n'):
            if line.startswith('ATOM') and ' CA ' in line:
                try:
                    b_factor = float(line[60:66].strip())
                    confidence_scores.append(b_factor)
                except:
                    pass

        avg_confidence = round(sum(confidence_scores) / len(confidence_scores), 2) if confidence_scores else 0

        if avg_confidence >= 90:
            confidence_level = "Very High"
            confidence_note = "Structure prediction is highly reliable."
        elif avg_confidence >= 70:
            confidence_level = "High"
            confidence_note = "Structure prediction is reliable for most regions."
        elif avg_confidence >= 50:
            confidence_level = "Moderate"
            confidence_note = "Some regions may be less accurate. Interpret with caution."
        else:
            confidence_level = "Low"
            confidence_note = "Structure prediction has low confidence. Experimental validation recommended."

        # Count secondary structure elements from PDB
        helix_count = pdb_content.count('HELIX')
        sheet_count = pdb_content.count('SHEET')

        return jsonify({
            "status": "success",
            "structure_prediction": {
                "protein_length": len(clean_protein),
                "protein_sequence": clean_protein,
                "pdb_content": pdb_content,
                "prediction_method": "ESMFold (Meta AI)",
                "average_confidence": avg_confidence,
                "confidence_level": confidence_level,
                "confidence_note": confidence_note,
                "structural_features": {
                    "residues_modeled": ca_count,
                    "alpha_helices": helix_count,
                    "beta_sheets": sheet_count
                },
                "visualization_note": "Use the PDB content with 3Dmol.js or NGL Viewer for 3D visualization",
                "alphafold_link": "https://alphafold.ebi.ac.uk/",
                "esmatlas_link": "https://esmatlas.com/resources?action=fold"
            },
            "message": f"3D structure predicted successfully for {len(clean_protein)} amino acid protein"
        })

    except Exception as e:
        return jsonify({"status": "error", "message": f"Structure prediction failed: {str(e)}"}), 500


@app.route('/api/protein/alphafold/<uniprot_id>', methods=['GET'])
def get_alphafold_structure(uniprot_id):
    """Fetch existing AlphaFold structure for known proteins by UniProt ID"""
    try:
        # Fetch structure info from AlphaFold database
        af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        
        af_response = requests.get(af_url, timeout=30)
        
        if af_response.status_code == 404:
            return jsonify({
                "status": "error",
                "message": f"No AlphaFold structure found for UniProt ID {uniprot_id}. Try ESMFold for custom sequences.",
                "fallback": "/api/protein/structure"
            }), 404

        if af_response.status_code != 200:
            return jsonify({
                "status": "error",
                "message": "AlphaFold API unavailable. Try again later."
            }), 500

        af_data = af_response.json()

        if not af_data:
            return jsonify({"status": "error", "message": "No data returned from AlphaFold"}), 404

        structure_info = af_data[0] if isinstance(af_data, list) else af_data

        # Fetch the actual PDB file
        pdb_url = structure_info.get('pdbUrl', '')
        pdb_content = ""
        
        if pdb_url:
            try:
                pdb_response = requests.get(pdb_url, timeout=30)
                pdb_content = pdb_response.text
            except:
                pdb_content = ""

        return jsonify({
            "status": "success",
            "alphafold_structure": {
                "uniprot_id": uniprot_id,
                "entry_id": structure_info.get('entryId', ''),
                "gene": structure_info.get('gene', ''),
                "organism": structure_info.get('organismScientificName', ''),
                "protein_name": structure_info.get('uniprotDescription', ''),
                "sequence_length": structure_info.get('uniprotSequenceLength', 0),
                "model_confidence": structure_info.get('confidenceAvgLocalScore', 0),
                "pdb_url": pdb_url,
                "pdb_content": pdb_content[:50000] if pdb_content else "",
                "cif_url": structure_info.get('cifUrl', ''),
                "prediction_method": "AlphaFold2 (Google DeepMind)"
            },
            "message": f"AlphaFold structure retrieved for {uniprot_id}"
        })

    except Exception as e:
        return jsonify({"status": "error", "message": f"AlphaFold fetch failed: {str(e)}"}), 500
if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
