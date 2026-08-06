window.MathJax = {
  tex: {
    // pymdownx.arithmatex wraps most $...$/$$...$$ spans server-side as \(...\)/\[...\],
    // but a handful of complex table-embedded equations arrive as raw HTML from the
    // manual's docx->markdown conversion and bypass that step. Recognizing the raw
    // dollar delimiters too lets MathJax pick those up on its own client-side scan.
    inlineMath: [["\\(", "\\)"], ["$", "$"]],
    displayMath: [["\\[", "\\]"], ["$$", "$$"]],
    processEscapes: true,
  },
};
