import openai
from transformers import pipeline
from Bio import Entrez
from pypdf import PdfReader

Entrez.email = "jaynaschick@gmail.com"

openai.api_key = "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"


def extract_text_from_pdf(pdf_path):
    """Extract text from a PDF file using pypdf."""
    with open(pdf_path, "rb") as file:
        reader = PdfReader(file)
        text = "\n".join([page.extract_text() for page in reader.pages if page.extract_text()])
    return text

def summarize_text(text):
    """Summarize medical text using Hugging Face transformers."""
    summarizer = pipeline("summarization", model="facebook/bart-large-cnn")
    return summarizer(text, max_length=150, min_length=50, do_sample=False)[0]['summary_text']

def search_pubmed(query):
    """Retrieve relevant research from PubMed."""
    Entrez.email = "your_email@example.com"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
    record = Entrez.read(handle)
    return [f"https://pubmed.ncbi.nlm.nih.gov/{id_}/" for id_ in record["IdList"]]

def explain_medical_terms(text):
    """Use OpenAI to explain medical terms in simple language."""
    openai.api_key = "your_openai_api_key"
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[{"role": "system", "content": "Explain medical terms in an easy-to-understand way."},
                  {"role": "user", "content": text}]
    )
    return response["choices"][0]["message"]["content"]

def main():
    choice = input("Enter '1' to upload a medical text file (PDF or TXT) or '2' to enter a diagnosis/treatment: ")
    
    if choice == '1':
        file_path = input("Enter the file path: ")
        if file_path.endswith(".pdf"):
            text = extract_text_from_pdf(file_path)
        elif file_path.endswith(".txt"):
            with open(file_path, "r", encoding="utf-8") as file:
                text = file.read()
        else:
            print("Unsupported file format.")
            return
    elif choice == '2':
        text = input("Enter the medical diagnosis or treatment: ")
    else:
        print("Invalid choice.")
        return
    
    print("\nSummarizing text...")
    summary = summarize_text(text)
    print("Summary:", summary)
    
    print("\nSearching for related research on PubMed...")
    pubmed_links = search_pubmed(text)
    print("Relevant research:")
    for link in pubmed_links:
        print(link)
    
    print("\nExplaining medical terms...")
    explanation = explain_medical_terms(text)
    print("Explanation:", explanation)

if __name__ == "__main__":
    main()
