"""Interface pequena: os calculos sao executados pelo programa em C."""
from pathlib import Path
import subprocess
import tkinter as tk
from tkinter import ttk, messagebox

PASTA = Path(__file__).resolve().parent
EXECUTAVEL = PASTA / "programa"
METODOS = {
    "Newton": "newton",
    "Ponto Fixo": "ponto-fixo",
    "Bisseção": "bissecao",
    "Gauss-Seidel": "gauss-seidel",
}

def executar():
    if not EXECUTAVEL.is_file():
        messagebox.showerror("Programa não encontrado", "Execute 'make' no terminal primeiro.")
        return
    try:
        resultado = subprocess.run(
            [str(EXECUTAVEL), METODOS[escolha.get()]],
            capture_output=True, text=True, timeout=10, check=False
        )
    except (OSError, subprocess.TimeoutExpired) as erro:
        messagebox.showerror("Erro", str(erro))
        return
    saida.configure(state="normal")
    saida.delete("1.0", "end")
    saida.insert("1.0", resultado.stdout or resultado.stderr)
    saida.configure(state="disabled")

janela = tk.Tk()
janela.title("Métodos Numéricos")
janela.geometry("660x420")
janela.minsize(530, 330)
conteudo = ttk.Frame(janela, padding=18)
conteudo.pack(fill="both", expand=True)

ttk.Label(conteudo, text="Métodos Numéricos", font=("Sans", 16, "bold")).pack(anchor="w", pady=(0, 12))
barra = ttk.Frame(conteudo)
barra.pack(fill="x", pady=(0, 12))
escolha = tk.StringVar(value="Gauss-Seidel")
ttk.Label(barra, text="Método:").pack(side="left", padx=(0, 10))
seletor = ttk.Combobox(barra, textvariable=escolha, values=list(METODOS), state="readonly", width=20)
seletor.pack(side="left", padx=(0, 10))
ttk.Button(barra, text="Executar", command=executar).pack(side="left")

saida = tk.Text(conteudo, font=("Monospace", 10), wrap="word", padx=12, pady=12,
                borderwidth=1, relief="solid", state="disabled")
saida.pack(fill="both", expand=True)
executar()
janela.mainloop()
