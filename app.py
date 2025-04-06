from flask import Flask, request, jsonify, render_template
from dotenv import load_dotenv
import os
from google import genai
import markdown  # Import the markdown library
from algo import ChemicalEquationBalancer

load_dotenv()

app = Flask(__name__)

balancer = ChemicalEquationBalancer()

def get_explanation(input_eq, balanced_eq):
    """
    Query Gemini's API to generate an explanation of how the equation was balanced.
    """
    prompt = (
        f"Explain in detail the steps involved in balancing the following chemical equation:\n\n"
        f"Input Equation: {input_eq}\n"
        f"Balanced Equation: {balanced_eq}\n\n"
        f"Explanation:"
    )
    
    client = genai.Client(api_key=os.getenv("GEMINI_API_KEY"))
    response = client.models.generate_content(
        model="gemini-2.0-flash",
        contents=prompt
    )
    
    # Return the explanation text
    return response.text.strip()

@app.route('/')
def home():
    # Render the HTML form
    return render_template('index.html')

@app.route('/balance', methods=['POST'])
def balance_form():
    # Get the equation from the submitted form
    equation = request.form.get('equation')
    try:
        balanced_equation = balancer.balance(equation)
        # Generate an explanation using Gemini's API
        explanation = get_explanation(equation, balanced_equation)
        # Convert markdown to HTML
        explanation_html = markdown.markdown(explanation)
        return render_template('result.html', equation=equation, balanced=balanced_equation, explanation=explanation_html)
    except Exception as e:
        return render_template('result.html', equation=equation, error=str(e))
    
if __name__ == '__main__':
    app.run(debug=True)