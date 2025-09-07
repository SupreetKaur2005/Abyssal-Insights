const carousel = document.querySelector('.carousel');
const cards = document.querySelectorAll('.species-card');
let index = 0;
const total = cards.length;

function showSlide(i) {
  index = (i + total) % total; 
  carousel.style.transform = `translateX(-${index * 100}%)`;
}

// Next / Prev buttons
document.querySelector('.next').addEventListener('click', () => showSlide(index + 1));
document.querySelector('.prev').addEventListener('click', () => showSlide(index - 1));

// Auto-slide every 10s
setInterval(() => {
  showSlide(index + 1);
}, 10000);
