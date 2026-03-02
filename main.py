from flask import Flask, request, jsonify
from flask_cors import CORS
from Bio import SeqIO
import io
import os
import requests

app = Flask(__name__)

CORS(app, origins="*", allow_headers=["Content-Type", "Authorization"], methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"])

def interpret_gc_content(gc_percentage, sequence_length, g_count, c_count, a_count, t_count):
    # GC Content Assessment
    if 40 <= gc_percentage <= 60:
        gc_status = "Normal"
        gc_interpretation = f"GC content of {gc_percentage}% falls within the normal range for human genomic sequences (40-60%). This indicates good sequence quality."
        gc_recommendation = "Sequence quality is acceptable. You may proceed with downstream analysis such as variant calling or alignment."
    elif gc_percentage < 40:
        gc_status = "AT-Rich"
        gc_interpretation = f"GC content of {gc_percentage}% is below the normal range (40-60%), indicating an AT-rich sequence. This may suggest a non-coding region, mitochondrial DNA variation, or low complexity region."
        gc_recommendation = "Flag for review. Consider verifying sample quality and re-sequencing if this result is unexpected for your target region."
    else:
        gc_status = "GC-Rich"
        gc_interpretation = f"GC content of {gc_percentage}% is above the normal range (40-60%), indicating a GC-rich sequence. This may suggest a promoter region, CpG island, or potential PCR amplification bias."
        gc_recommendation = "Flag for review. Verify that PCR conditions were optimized for GC-rich templates. Consider checking for amplification bias."

    # Sequence Length Assessment
    if sequence_length < 100:
        length_status = "Very Short"
        length_interpretation = "Sequence length is very short (under 100 bp). This may indicate an incomplete or degraded sample."
        length_recommendation = "Re-sequence the sample to obtain a higher coverage and longer read length."
    elif sequence_length < 1000:
        length_status = "Short"
        length_interpretation = f"Sequence length of {sequence_length} bp is relatively short. Suitable for targeted sequencing analysis but limited for whole genome studies."
        length_recommendation = "Suitable for targeted gene panels. For whole genome analysis, consider longer read sequencing."
    elif sequence_length < 50000:
        length_status = "Standard"
        length_interpretation = f"Sequence length of {sequence_length} bp is within standard range for genomic analysis. This is consistent with mitochondrial genome or targeted sequencing data."
        length_recommendation = "Proceed with standard bioinformatics pipeline for alignment and variant calling."
    else:
        length_status = "Long"
        length_interpretation = f"Sequence length of {sequence_length} bp indicates a large genomic region or whole chromosome data. High coverage and computational resources recommended."
        length_recommendation = "Use high-performance computing pipeline for alignment. Consider splitting into smaller chunks for faster processing."

    # Nucleotide Balance Assessment
    total = g_count + c_count + a_count + t_count
    a_pct = round((a_count / total) * 100, 2) if total > 0 else 0
    t_pct = round((t_count / total) * 100, 2) if total > 0 else 0
    g_pct = round((g_count / total) * 100, 2) if total > 0 else 0
    c_pct = round((c_count / total) * 100, 2) if total > 0 else 0

    # GC Skew
    gc_skew = round((g_count - c_count) / (g_count + c_count), 4) if (g_count + c_count) > 0 else 0
    if gc_skew > 0.1:
        skew_interpretation = "Positive GC skew detected, suggesting leading strand bias. Common in bacterial genomes and mitochondrial DNA."
    elif gc_skew < -0.1:
        skew_interpretation = "Negative GC skew detected, suggesting lagging strand bias."
    else:
        skew_interpretation = "GC skew is balanced, indicating no significant strand bias."

    # Overall Quality Score (0-100)
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
        overall_summary = "This sample demonstrates good quality metrics with minor flags. Results are reliable for most analyses."
    elif quality_score >= 40:
        overall_quality = "Fair"
        overall_summary = "This sample has some quality concerns. Review flagged metrics before proceeding with clinical interpretation."
    else:
        overall_quality = "Poor"
        overall_summary = "This sample has significant quality issues. Re-sequencing is strongly recommended before clinical use."

    return {
        "gc_status": gc_status,
        "gc_interpretation": gc_interpretation,
        "gc_recommendation": gc_recommendation,
        "length_status": length_status,
        "length_interpretation": length_interpretation,
        "length_recommendation": length_recommendation,
        "nucleotide_percentages": {
            "A": a_pct,
            "T": t_pct,
            "G": g_pct,
            "C": c_pct
        },
        "gc_skew": gc_skew,
        "skew_interpretation": skew_interpretation,
        "quality_score": quality_score,
        "overall_quality": overall_quality,
        "overall_summary": overall_summary
    }


@app.route('/', methods=['GET'])
def home():
    return jsonify({
        "app_name": "Genomirates",
        "status": "online",
        "message": "Genomirates Bioinformatics Engine is running.",
        "version": "2.0",
        "endpoints": {
            "analyze": "POST /api/calc — Analyze a genomic sequence",
        }
    })


@app.route('/api/calc', methods=['GET', 'POST'])
def analyze_sequence():
    if request.method == 'GET':
        return jsonify({
            "app_name": "Genomirates",
            "status": "online",
            "message": "Send a POST request with a 'sequence' or 'file_url' field to analyze DNA."
        })

    try:
        data = request.get_json()
        if not data:
            return jsonify({
                "app_name": "Genomirates",
                "status": "error",
                "message": "Missing JSON body"
            }), 400

        fasta_text = ""
        if 'file_url' in data:
            try:
                response = requests.get(data['file_url'], timeout=30)
                response.raise_for_status()
                fasta_text = response.text
            except Exception as e:
                return jsonify({
                    "app_name": "Genomirates",
                    "status": "error",
                    "message": f"Failed to download file: {str(e)}"
                }), 400
        elif 'sequence' in data:
            fasta_text = data['sequence']
        else:
            return jsonify({
                "app_name": "Genomirates",
                "status": "error",
                "message": "Missing 'sequence' or 'file_url' field in JSON body"
            }), 400

        fasta_io = io.StringIO(fasta_text)
        try:
            record = next(SeqIO.parse(fasta_io, "fasta"))
        except StopIteration:
            return jsonify({
                "app_name": "Genomirates",
                "status": "error",
                "message": "Invalid FASTA format"
            }), 400

        sequence_str = str(record.seq).upper()
        if len(sequence_str) == 0:
            return jsonify({
                "app_name": "Genomirates",
                "status": "error",
                "message": "Empty sequence detected"
            }), 400

        sequence_length = len(sequence_str)
        g_count = sequence_str.count('G')
        c_count = sequence_str.count('C')
        a_count = sequence_str.count('A')
        t_count = sequence_str.count('T')
        gc_content = round(((g_count + c_count) / sequence_length) * 100, 2) if sequence_length > 0 else 0

        interpretation = interpret_gc_content(gc_content, sequence_length, g_count, c_count, a_count, t_count)

        return jsonify({
            "app_name": "Genomirates",
            "status": "success",
            "analysis_results": {
                "sample_id": record.id,
                "sequence_length": sequence_length,
                "gc_content_percentage": gc_content,
                "g_count": g_count,
                "c_count": c_count,
                "a_count": a_count,
                "t_count": t_count,
                "nucleotide_percentages": interpretation["nucleotide_percentages"],
                "gc_skew": interpretation["gc_skew"],
                "gc_status": interpretation["gc_status"],
                "gc_interpretation": interpretation["gc_interpretation"],
                "gc_recommendation": interpretation["gc_recommendation"],
                "length_status": interpretation["length_status"],
                "length_interpretation": interpretation["length_interpretation"],
                "length_recommendation": interpretation["length_recommendation"],
                "skew_interpretation": interpretation["skew_interpretation"],
                "quality_score": interpretation["quality_score"],
                "overall_quality": interpretation["overall_quality"],
                "overall_summary": interpretation["overall_summary"]
            },
            "message": f"Successfully analyzed {record.id}"
        })

    except Exception as e:
        return jsonify({
            "app_name": "Genomirates",
            "status": "error",
            "message": f"Processing failed: {str(e)}"
        }), 500


if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
    
