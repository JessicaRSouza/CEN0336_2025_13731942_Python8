#!/usr/bin/env python3

#Importando o módulo "re" para encontrar futuramente as ocorrências de metionina
import re

#Importando a tabela de tradução de aminoácidos
tabela_de_traducao = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}


#Solicitando o nome do arquivo ao usuário
arquivo = input("Digite o nome do arquivo FASTA: ")

#Criando os dicionários (vazios por agora)
fastaDict = {}  #Armazena os identificadores e suas sequências

#Abrindo o arquivo FASTA (leitura) e os quatro arquivos de saída (escrita)
with open(arquivo, "rt") as fasta, open("Python_08.codons-6frames.nt", "w") as arquivo_saida, open("Python_08.translated.aa", "w") as arquivo_aa, open("Python_08.translated-longest.aa", "w") as arq_pept_longo, open("Python_08.orf-longest.nt", "w") as arq_orf:
	for line in fasta:
		line = line.rstrip()  #Remove quebras de linha e espaços ao fim

		if line.startswith(">"):  #Identifica o cabeçalho de cada sequência
			identificador = line.split(" ")[0].replace(">","")  #Extrai o identificador da sequência, dividindo a linha pelo primeiro espaço. Também remove o caractere ">", para uma saída mais limpa.
			if identificador not in fastaDict.keys():
				#Caso o identificador ainda não seja uma chave no dicionário, inicia o valor com uma string vazia para acumular a sequência
				fastaDict[identificador] = ""

		else:
			#Adiciona o fragmento de sequência à sequência total correspondente ao identificador
			#Certifica que a sequência está em maiúsculas
			fastaDict[identificador] += line.upper()


	#Nesse ponto, o arquivo FASTA foi lido e o dicionário fastaDict está completo
	#Agora, iteramos pelo dicionário para processar cada sequência e escrever no arquivo de saída
	for identificador, sequencia in fastaDict.items():

		#Criando as variáveis que guardarão o peptídeo mais longo dos 6 quadros
		pept_mais_longo = ""  #Guarda o peptídeo mais longo
		codons_pept_mais_longo = ""  #Guarda os códons do peptídeo mais longo

		#Obtendo a sequência complementar reversa
		seq_reversa = sequencia[::-1]  #Obtém o reverso ao inverter a string
		seq_comp_rev = seq_reversa.replace("T", "a").replace("A", "t").replace("G", "c").replace("C", "g").upper()  #Obtém o complemento a partir da utilização de minúsculas temporárias

		#Criando uma lista das sequências a serem processadas -> cria-se uma "lista de tarefas"
		#Cada  tarefa é uma tupla (par) contendo: (a_sequência_para_processar, o_número_base_do_quadro)
		seqs_processar = [
			(sequencia, 0),  #Tarefa 1: Processar a "sequência original", usando 0 como "base" do quadro
			(seq_comp_rev, 3)  # Tarefa 1: Processar a "sequência complementar reversa", usando 3 como "base" do quadro
		]

		#Criando um loop para iterar pelas sequências a processar
		for seq_atual, base_quadro in seqs_processar:

			#Criando um novo loop para iterar pelos distintos quadros de leitura
			for quadro_local in range (1,4):

				inicio_index = quadro_local-1  #Calcula o índice inicial (0, 1 ou 2) para o "slice" da string

				quadro_global = quadro_local + base_quadro  #Calcula o quadro "global" somando o quadro local (1,2,3) ao número 'base' (0 ou 3)

				#Escrevendo o cabeçalho formatado para ambos os arquivos
				arquivo_saida.write(f"{identificador}-frame-{quadro_global}-codons\n")
				arquivo_aa.write(f"{identificador}-frame-{quadro_global}-translation\n")

				codons = []  #Cria uma lista vazia para armazenar os códons da sequência atual
				aminoacidos = []  #Cria uma segunda lista vazia para armazenar os aminoácidos

				#Iterando pela string da sequência, começando do "inicio_index" e pulando de 3 em 3 bases (quadro de leitura 1)
				for i in range(inicio_index, len(seq_atual), 3):
					codon = seq_atual[i:i+3]  #Extrai o slice de 3 bases
					if len(codon) == 3:  #Garante que o códon esteja completo (ignora "restos" no final da sequência)
						codons.append(codon)  #Adiciona o códon completo à lista

						#Traduzindo o códon para aminoácido
						aa = tabela_de_traducao.get(codon, "X")
						aminoacidos.append(aa)  #Adiciona o aminoácido à lista

				#Juntando todos os códons da lista, separados por um espaço
				linha_codons_formatada = " ".join(codons)

				#Escrevendo a linha de códons formatada no arquivo 1 e adicionando uma linha em branco para separar as entradas no arquivo de saída
				arquivo_saida.write(f"{linha_codons_formatada}\n")
				arquivo_saida.write("\n")

				#Juntando todos os aminoácidos da lista
				linha_aa_formatada = "".join(aminoacidos)

				#Escrevendo a linha de aminoácidos formatada no arquivo 2 e adicionando uma linha em branco para separar as entradas no arquivo de saída
				arquivo_aa.write(f"{linha_aa_formatada}\n")
				arquivo_aa.write("\n")

				#Procurando por todas as metioninas ("M"), que são os inícios de peptídeos
				for match in re.finditer("M", linha_aa_formatada):

					m_index = match.start()  #Determina o índice onde o "M" está
					sub_seq_comeco = linha_aa_formatada[m_index:]  #Fatia da string que começa a partir desse "M"
					sub_seq_fim = sub_seq_comeco.find("*")  #Procura o primeiro stop codon ("*") dentro dessa fatia

					#Criando variáveis temporárias para o peptídeo e códons atuais
					pept_atual = ""
					codons_atuais = []

					if sub_seq_fim == -1:
						#Determina o que fazer quando não encontra stop codon
						pept_atual = sub_seq_comeco  #O peptídeo segue até o fim da sequência
						codons_atuais = codons[m_index:]  #Pega os códons correspondentes (do "M" até o fim)

					else:
						#Determina o que fazer quando o stop codon é encontrado
						pept_atual = sub_seq_comeco[:sub_seq_fim]  #O peptídeo vai do "M" até o stop codon
						fim_index = m_index + sub_seq_fim  #Calcula o índice final para fatiar a lista de códons
						codons_atuais = codons[m_index:fim_index]  #Pega os códons correspondentes (do "M" até o stop codon)

					#Verificando se o peptídeo atual é o mais longo
					if len(pept_atual) > len(pept_mais_longo):
						pept_mais_longo = pept_atual  #Caso seja, ele é definido como o novo peptídeo mais longo
						codons_pept_mais_longo = "".join(codons_atuais)  #Define os códons do peptídeo mais longo

		#Após checar os 6 quadros, escreve os maiores nos arquivos 3 e 4

		#Escrevendo o peptídeo mais longo no arquivo 3
		arq_pept_longo.write(f">{identificador}-longest-peptide-M-to-Stop\n")
		arq_pept_longo.write(f"{pept_mais_longo}\n\n")
		#Escrevendo os códons do peptídeo mais longo no arquivo 4
		arq_orf.write(f">{identificador}-longest-orf-nucleotides\n")
		arq_orf.write(f"{codons_pept_mais_longo}\n\n")

#Emitindo uma mensagem final para o usuário, indicando que o processo terminou
print("Processamento concluído. Verifique os arquivos \"Python_08.codons-6frames.nt\", \"Python_08.translated.aa\" e \"Python_08.translated-longest.aa\".")
