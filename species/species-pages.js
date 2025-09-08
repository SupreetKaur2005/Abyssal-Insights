const facts = document.querySelectorAll('.fact');
let factIndex = 0;

function showFact(i) {
  facts.forEach((fact, idx) => {
    fact.style.transform = `translateX(${100 * (idx - i)}%)`;
  });
}

setInterval(() => {
  factIndex = (factIndex + 1) % facts.length;
  showFact(factIndex);
}, 5000);

showFact(factIndex);
